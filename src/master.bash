#!/bin/bash

# This is the master script to run all pieces of the pipeline to parse clinvar into a tab-delimited file
# example usage: bsub -q priority -o cv.o -e cv.e -J clinvar -R rusage[mem=32] "cd $workdir; . ./private_paths.bash; cd clinvar; ./master.bash"
# output: clinvar.tsv

# required environment variables:
# $vt - path to vt
# $b37ref - path to a b37 .fa file

# download latest clinvar XML and tab-delimited summary
rm ClinVarFullRelease_00-latest.xml.gz
rm variant_summary.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# extract the GRCh37 coordinates, mutant allele, MeasureSet ID and PubMed IDs from it
./parse_clinvar_xml.py -x ClinVarFullRelease_00-latest.xml.gz -o clinvar_table_raw.tsv
# the above takes about 3 hours (with submitters + traits, old version was ~1h without), so recommend submitting as a job:
# bsub -q priority -R rusage[mem=64] -oo cvxml.o -eo cvxml.e -J cvxml "./parse_clinvar_xml.py -x ClinVarFullRelease_00-latest.xml.gz -o clinvar_table_raw.tsv"

# sort the table
cat clinvar_table_raw.tsv | head -1 > clinvar_table_sorted.tsv # header row
cat clinvar_table_raw.tsv | tail -n +2 | egrep -v "^[XYM]" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv # numerically sort chroms 1-22
cat clinvar_table_raw.tsv | tail -n +2 | egrep "^[XYM]" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv # lexicographically sort non-numerical chroms at end

# de-duplicate records
./dedup_clinvar.py < clinvar_table_sorted.tsv > clinvar_table_dedup.tsv

# remove hyphens and add one base pair of sequence context for indels to make them vt normalize-compatible. runs in <1 min
./add_context_to_indels.py -R $b37ref < clinvar_table_dedup.tsv > clinvar_table_dedup_context.tsv

# create a VCF and run vt-normalize
echo "##fileformat=VCFv4.1" > to_normalize.vcf
cat clinvar_table_dedup_context.tsv | awk -v FS="\t" -v OFS="\t" 'BEGIN {print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}; NR>1 {print $1,$2,".",$3,$4,".",".","."}' >> to_normalize.vcf
$vt normalize to_normalize.vcf -r $b37ref -o normalized.vcf

# paste the newly normalized CHROM POS REF ALT and the additional fields back together again
echo -e "chrom\tpos\tref\talt\tmut\tmeasureset_id\tall_submitters\tall_traits\tall_pmids" > clinvar_table_dedup_normalized.tsv
mkfifo coordinates
cat normalized.vcf | tail -n +5 | cut -f1,2,4,5 > coordinates &
cat clinvar_table_dedup_context.tsv | tail -n +2 | cut -f5-9 | paste coordinates - >> clinvar_table_dedup_normalized.tsv

# join information from the tab-delimited summary to the normalized genomic coordinates
Rscript join_data.R

# now sort again by genomic coordinates (because R's merge function ruins this)
cat clinvar_combined.tsv | head -1 > clinvar_combined_sorted.tsv # header row
cat clinvar_combined.tsv | tail -n +2 | egrep -v "^[XYM]" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> clinvar_combined_sorted.tsv # numerically sort chroms 1-22
cat clinvar_combined.tsv | tail -n +2 | egrep "^[XYM]" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_combined_sorted.tsv # lexicographically sort non-numerical chroms at end

# now de-dup _again_, because the tab-delimited summary contains dups
./dedup_clinvar.py < clinvar_combined_sorted.tsv > clinvar_combined_sorted_dedup.tsv

# create a text file
cp clinvar_combined_sorted_dedup.tsv clinvar.tsv
# placeholder in case it turns out we need to add more steps

# to do: create a VCF
# echo "##fileformat=VCFv4.1" > clinvar.vcf
# cat clinvar.tsv | awk -v FS="\t" -v OFS="\t" 'BEGIN {print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}; NR==1 do something {NR>1 {print $1,$2,".",$3,$4,".",".","MUT="$5";MEASURESET_ID="$6";PMIDS="$7}' >> to_normalize.vcf

# clean up
# rm clinvar_table_raw.tsv # saving this for now b/c it takes an hour to generate and i am not sure this script works yet
rm clinvar_table_sorted.tsv
rm clinvar_table_dedup.tsv
rm clinvar_table_dedup_context.tsv
rm to_normalize.vcf
rm normalized.vcf
rm clinvar_table_dedup_normalized.tsv
rm clinvar_combined_sorted.tsv

