"""
Script for generating a new clinvar table with the gnomAD fields below added to each clinvar variant that's in gnomAD exomes or genomes
"""
import argparse
from collections import defaultdict
import gzip
import pysam
import sys

NEEDED_GNOMAD_FIELDS = [ 'Filter',  # whether the variant is PASS
 'AC', 'AN', 'AF', 'DP','Hom',
 'AC_AFR', 'AC_AMR', 'AC_ASJ', 'AC_EAS', 'AC_SAS', 'AC_FIN', 'AC_NFE', 'AC_OTH', 
 'AN_AFR', 'AN_AMR', 'AN_ASJ', 'AN_EAS', 'AN_SAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 
 'AC_AFR', 'AF_AMR', 'AF_ASJ', 'AF_EAS', 'AF_SAS', 'AF_FIN', 'AF_NFE', 'AF_OTH', 
 'AC_Male', 'AC_Female', 'AN_Male', 'AN_Female',
 'Hom_AFR', 'Hom_AMR', 'Hom_ASJ', 'Hom_EAS', 'Hom_SAS', 'Hom_FIN', 'Hom_NFE', 'Hom_OTH', 
 'Hemi_AFR', 'Hemi_AMR', 'Hemi_ASJ', 'Hemi_EAS', 'Hemi_SAS', 'Hemi_FIN', 'Hemi_NFE', 'Hemi_OTH', 
 'Hom_Male', 'Hom_Female',
 'AS_RF', 'AS_FilterStatus',
 'AC_POPMAX', 'AN_POPMAX', 'AF_POPMAX', 'POPMAX', 
]

NEEDED_GNOMAD_FIELDS_SET = set(NEEDED_GNOMAD_FIELDS)
GNOMAD_EMPTY_COLUMN_VALUES = ['']*len(NEEDED_GNOMAD_FIELDS_SET)

p = argparse.ArgumentParser()
p.add_argument("-i", "--clinvar-table", help="Clinvar .tsv", required=True)
g = p.add_mutually_exclusive_group(required=True)
g.add_argument("-ge", "--gnomad-exomes-vcf", dest="gnomad_sites_vcf", help="gnomAD exomes VCF directory")
g.add_argument("-gg", "--gnomad-genomes-vcf", dest="gnomad_sites_vcf", help="gnomAD genomes VCF directory")
args = p.parse_args()

counts = defaultdict(int)

def get_gnomad_column_values(gnomad_f, chrom, pos, ref, alt):
    """Retrieves the gnomAD vcf row corresponding to the given chrom, pos, ref, alt, and extracts the column values listed in NEEDED_GNOMAD_FIELDS

    Args:
      gnomad_f: A pysam.TabixFile object corresponding to the gnomAD exomes or genomes vcf
      chrom: chromosome (eg. '1')
      pos: the minrepped clinvar variant position
      ref: the minrepped clinvar ref allele
      ref: the minrepped clinvar alt allele

    Return:
      A list of chrom, pos, ref, alt
    """

    if chrom == 'MT':
        return GNOMAD_EMPTY_COLUMN_VALUES

    counts['total_clinvar_variants'] += 1

    # retrieve gnomAD variant - pysam.fetch(..) sometimes returns more than 1 vcf record, so need to filter here
    position_found = False
    gnomad_alt_alleles = []
    for gnomad_vcf_row in gnomad_f.fetch(chrom, pos-1, pos):
        gnomad_row_fields = gnomad_vcf_row.split('\t')
        if str(pos) !=  gnomad_row_fields[1]:
            continue
        position_found = True
        gnomad_ref_allele = gnomad_row_fields[3]
        gnomad_alt_allele = gnomad_row_fields[4]
        if "," in gnomad_alt_allele:
            raise Exception("Found multiallelic variant: %s. Expecting an gnomAD VCF that has been decomposed / normalized with vt." % "-".join(gnomad_vcf_row_fields[0:5]))

        if ref == gnomad_ref_allele and alt == gnomad_alt_allele:
            counts['clinvar_variants_with_matching_position_and_matching_allele'] += 1
            break
        gnomad_alt_alleles.append(gnomad_alt_allele)
    else:
        if not position_found:
            counts['clinvar_variants_with_no_matching_position_in_gnomad'] += 1
        else:
            if len(ref) + len(alt) + len(gnomad_ref_allele) + len(gnomad_alt_allele) > 4:
                counts['clinvar_indel_with_no_matching_allele_in_gnomad'] += 1                
            elif ref != gnomad_ref_allele and alt != gnomad_alt_allele:
                counts['clinvar_snp_with_mismatching_ref_and_alt_allele_in_gnomad'] += 1   
            elif ref != gnomad_ref_allele:
                counts['clinvar_snp_with_mismatching_ref_allele_in_gnomad'] += 1
            elif alt != gnomad_alt_allele:
                counts['clinvar_snp_with_mismatching_alt_allele_in_gnomad'] += 1
            else:
                counts['clinvar_snp_with_unknown_mismatch'] += 1

            sys.stderr.write("WARNING: gnomAD variant %s:%s (http://gnomad.broadinstitute.org/variant/%s-%s-%s-%s) - gnomAD alleles (%s:%s %s>%s) mismatch the clinvar allele (%s:%s %s>%s)\n" % (chrom, pos, chrom, pos, gnomad_row_fields[3], gnomad_row_fields[4], chrom, pos, gnomad_ref_allele, ",".join(gnomad_alt_alleles), chrom, pos, ref, alt))

        return GNOMAD_EMPTY_COLUMN_VALUES

    filter_value = gnomad_row_fields[6]
    info_fields = [('Filter', filter_value)] + [tuple(kv.split('=')) for kv in gnomad_row_fields[7].split(';')]
    info_fields = filter(lambda kv: kv[0] in NEEDED_GNOMAD_FIELDS_SET, info_fields)
    info_fields = dict(info_fields)
    gnomad_column_values = [info_fields.get(k, '') for k in NEEDED_GNOMAD_FIELDS]

    # check that the clinvar alt allele matches (one of the) gnomAD alt allele(s)    
    #if len(alt_alleles) > 1:
    #    # select the AC/AN numbers corresponding to the specific alt allele
    #    alt_allele_index = alt_alleles.index(alt)    
    #    gnomad_column_values = map(lambda x: x.split(",")[alt_allele_index] if "," in x else x, gnomad_column_values)

    return gnomad_column_values


gnomad_f = pysam.TabixFile(args.gnomad_sites_vcf)
clinvar_f = gzip.open(args.clinvar_table) if args.clinvar_table.endswith('.gz') else open(args.clinvar_table)
clinvar_header = next(clinvar_f).rstrip('\n').split('\t')
clinvar_with_gnomad_header = clinvar_header + NEEDED_GNOMAD_FIELDS
print("\t".join(clinvar_with_gnomad_header))
for i, clinvar_row in enumerate(clinvar_f):
    clinvar_fields = clinvar_row.rstrip('\n').split('\t')
    clinvar_dict = dict(zip(clinvar_header, clinvar_fields))

    chrom = clinvar_dict['chrom']
    pos = int(clinvar_dict['pos'])
    ref = clinvar_dict['ref']
    alt = clinvar_dict['alt']
    gnomad_column_values = get_gnomad_column_values(gnomad_f, chrom, pos, ref, alt)
    
    print("\t".join(clinvar_fields + gnomad_column_values))

for k, v in counts.items():
    sys.stderr.write("%30s: %s\n" % (k, v))
