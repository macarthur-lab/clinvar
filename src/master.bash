#!/bin/bash

# This is the master script to run all pieces of the pipeline to parse clinvar into a tab-delimited file
# example usage: bsub -q priority -o cv.o -e cv.e -J clinvar -R rusage[mem=32] "cd $workdir; . ./private_paths.bash; cd clinvar; ./master.bash"
# output: clinvar.tsv

# required environment variables:
# $b37ref - path to a b37 .fa file

if [ -z "$b37ref" ]; then
   echo "\$b37ref must be set to a b37 .fa file"
   exit -1;
fi;

set -x #echo on

# download latest clinvar XML and tab-delimited summary
rm ClinVarFullRelease_00-latest.xml.gz
rm variant_summary.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# extract the GRCh37 coordinates, mutant allele, MeasureSet ID and PubMed IDs from it
python parse_clinvar_xml.py -x ClinVarFullRelease_00-latest.xml.gz -o clinvar_table_raw.tsv
# the above now takes about 20 minutes. if you want to submit as a job instead:
# bsub -q priority -R rusage[mem=8] -oo cvxml.o -eo cvxml.e -J cvxml "./parse_clinvar_xml.py -x ClinVarFullRelease_00-latest.xml.gz -o clinvar_table_raw.tsv"

# sort the table
rm clinvar_table_sorted.tsv
cat clinvar_table_raw.tsv | head -1 > clinvar_table_sorted.tsv # header row
cat clinvar_table_raw.tsv | tail -n +2 | egrep -v "^[XYM]" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv # numerically sort chroms 1-22
cat clinvar_table_raw.tsv | tail -n +2 | egrep "^[XYM]" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv # lexicographically sort non-numerical chroms at end

# de-duplicate records
rm clinvar_table_dedup.tsv
rm clinvar_table_dedup_context.tsv
python dedup_clinvar.py < clinvar_table_sorted.tsv > clinvar_table_dedup.tsv

# normalize (convert to minimal representation and left-align)
# the normalization code is in a different repo (useful for more than just clinvar) so here I just wget it:
rm to_normalize.vcf
rm normalized.vcf
rm clinvar_table_dedup_normalized.tsv
wget -N https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py
python normalize.py -R $b37ref < clinvar_table_dedup.tsv > clinvar_table_dedup_normalized.tsv

# join information from the tab-delimited summary to the normalized genomic coordinates
Rscript join_data.R

# now sort again by genomic coordinates (because R's merge function ruins this)
rm clinvar_combined_sorted.tsv
cat clinvar_combined.tsv | head -1 > clinvar_combined_sorted.tsv # header row
cat clinvar_combined.tsv | tail -n +2 | egrep -v "^[XYM]" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> clinvar_combined_sorted.tsv # numerically sort chroms 1-22
cat clinvar_combined.tsv | tail -n +2 | egrep "^[XYM]" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_combined_sorted.tsv # lexicographically sort non-numerical chroms at end

# now de-dup _again_, because the tab-delimited summary contains dups
python dedup_clinvar.py < clinvar_combined_sorted.tsv > clinvar_combined_sorted_dedup.tsv

# create a text file
cp clinvar_combined_sorted_dedup.tsv clinvar.tsv
bgzip -c clinvar.tsv > clinvar.tsv.gz  # create compressed version
tabix -S 1 -s 1 -b 2 -e 2 clinvar.tsv.gz


# create vcf
python clinvar_table_to_vcf.py -o clinvar.vcf clinvar.tsv
bgzip -c clinvar.vcf > clinvar.vcf.gz  # create compressed version
tabix clinvar.vcf.gz


# clean up
rm clinvar_table_raw.tsv
rm clinvar_table_sorted.tsv
rm clinvar_table_dedup.tsv
rm clinvar_table_dedup_context.tsv
rm to_normalize.vcf
rm normalized.vcf
rm clinvar_table_dedup_normalized.tsv
rm clinvar_combined_sorted.tsv