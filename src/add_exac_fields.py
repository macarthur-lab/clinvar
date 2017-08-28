"""
Script for generating a new clinvar table with the ExAC fields below added to each clinvar variant that's in ExAC
"""
import argparse
from collections import defaultdict
import gzip
import pysam
import sys

NEEDED_EXAC_FIELDS = [ 'Filter',  # whether the variant is PASS
 'AC', 'AC_Het', 'AC_Hom', 'AC_Adj', 'AN', 'AN_Adj', 'AF', 
 'AC_AFR', 'AC_AMR', 'AC_EAS', 'AC_FIN', 'AC_NFE', 'AC_OTH', 'AC_SAS', 
 'AN_AFR', 'AN_AMR', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AN_SAS',
 'Het_AFR', 'Het_AMR', 'Het_EAS', 'Het_FIN', 'Het_NFE', 'Het_OTH', 'Het_SAS', 
 'Hom_AFR', 'Hom_AMR', 'Hom_EAS', 'Hom_FIN', 'Hom_NFE', 'Hom_OTH', 'Hom_SAS', 
 'AC_MALE', 'AC_FEMALE', 'AN_MALE', 'AN_FEMALE', 'AC_CONSANGUINEOUS', 'AN_CONSANGUINEOUS', 'Hom_CONSANGUINEOUS', 
 'AC_POPMAX', 'AN_POPMAX', 'POPMAX']

NEEDED_EXAC_FIELDS_SET = set(NEEDED_EXAC_FIELDS)
EXAC_EMPTY_COLUMN_VALUES = ['']*len(NEEDED_EXAC_FIELDS_SET)

p = argparse.ArgumentParser()
p.add_argument("-i", "--clinvar-table", help="Clinvar .tsv", required=True)
p.add_argument("-e", "--exac-sites-vcf", help="ExAC sites VCF", required=True)
args = p.parse_args()

counts = defaultdict(int)

def get_exac_column_values(exac_f, chrom, pos, ref, alt):
    """Retrieves the ExAC vcf row corresponding to the given chrom, pos, ref, alt, and extracts the column values listed in NEEDED_EXAC_FIELDS

    Args:
      exac_f: A pysam.TabixFile object corresponding to the ExAC vcf
      chrom: chromosome (eg. '1')
      pos: the minrepped clinvar variant position
      ref: the minrepped clinvar ref allele
      ref: the minrepped clinvar alt allele

    Return:
      A list of chrom, pos, ref, alt
    """

    if chrom == 'MT':
        return EXAC_EMPTY_COLUMN_VALUES

    counts['total_clinvar_variants'] += 1

    # retrieve ExAC variant - pysam.fetch(..) sometimes returns more than 1 vcf record, so need to filter here
    position_found = False
    exac_alt_alleles = []
    for exac_vcf_row in exac_f.fetch(chrom, pos-1, pos):
        exac_row_fields = exac_vcf_row.split('\t')
        if str(pos) !=  exac_row_fields[1]:
            continue
        position_found = True
        exac_ref_allele = exac_row_fields[3]
        exac_alt_allele = exac_row_fields[4]
        if "," in exac_alt_allele:
            raise Exception("Found multiallelic variant: %s. Expecting an ExAC VCF that has been decomposed / normalized with vt." % "-".join(exac_vcf_row_fields[0:5]))

        if ref == exac_ref_allele and alt == exac_alt_allele:
            counts['clinvar_variants_with_matching_position_and_matching_allele'] += 1
            break
        exac_alt_alleles.append(exac_alt_allele)
    else:
        if not position_found:
            counts['clinvar_variants_with_no_matching_position_in_exac'] += 1
        else:
            if len(ref) + len(alt) + len(exac_ref_allele) + len(exac_alt_allele) > 4:
                counts['clinvar_indel_with_no_matching_allele_in_exac'] += 1                
            elif ref != exac_ref_allele and alt != exac_alt_allele:
                counts['clinvar_snp_with_mismatching_ref_and_alt_allele_in_exac'] += 1   
            elif ref != exac_ref_allele:
                counts['clinvar_snp_with_mismatching_ref_allele_in_exac'] += 1
            elif alt != exac_alt_allele:
                counts['clinvar_snp_with_mismatching_alt_allele_in_exac'] += 1
            else:
                counts['clinvar_snp_with_unknown_mismatch'] += 1

            sys.stderr.write("WARNING: ExAC variant %s:%s (http://exac.broadinstitute.org/variant/%s-%s-%s-%s) - ExAC alleles (%s:%s %s>%s) mismatch the clinvar allele (%s:%s %s>%s)\n" % (chrom, pos, chrom, pos, exac_row_fields[3], exac_row_fields[4], chrom, pos, exac_ref_allele, ",".join(exac_alt_alleles), chrom, pos, ref, alt))

        return EXAC_EMPTY_COLUMN_VALUES

    filter_value = exac_row_fields[6]
    info_fields = [('Filter', filter_value)] + [tuple(kv.split('=')) for kv in exac_row_fields[7].split(';')]
    info_fields = filter(lambda kv: kv[0] in NEEDED_EXAC_FIELDS_SET, info_fields)
    info_fields = dict(info_fields)
    try:
        exac_column_values = [info_fields[k] for k in NEEDED_EXAC_FIELDS]
    except Exception, e:
        print([('Filter', filter_value)] + [tuple(kv.split('=')) for kv in exac_row_fields[7].split(';')])
        print(info_fields)
        raise ValueError("ERROR: unable to parse INFO fields in row: %s.  %s" % (exac_row_fields, e)) 

    # check that the clinvar alt allele matches (one of the) ExAC alt allele(s)    
    #if len(alt_alleles) > 1:
    #    # select the AC/AN numbers corresponding to the specific alt allele
    #    alt_allele_index = alt_alleles.index(alt)    
    #    exac_column_values = map(lambda x: x.split(",")[alt_allele_index] if "," in x else x, exac_column_values)

    return exac_column_values


exac_f = pysam.TabixFile(args.exac_sites_vcf)
clinvar_f = gzip.open(args.clinvar_table) if args.clinvar_table.endswith('.gz') else open(args.clinvar_table)
clinvar_header = next(clinvar_f).rstrip('\n').split('\t')
clinvar_with_exac_header = clinvar_header + NEEDED_EXAC_FIELDS
print("\t".join(clinvar_with_exac_header))
for i, clinvar_row in enumerate(clinvar_f):
    clinvar_fields = clinvar_row.rstrip('\n').split('\t')
    clinvar_dict = dict(zip(clinvar_header, clinvar_fields))

    chrom = clinvar_dict['chrom']
    pos = int(clinvar_dict['pos'])
    ref = clinvar_dict['ref']
    alt = clinvar_dict['alt']
    exac_column_values = get_exac_column_values(exac_f, chrom, pos, ref, alt)
    
    print("\t".join(clinvar_fields + exac_column_values))

for k, v in counts.items():
    sys.stderr.write("%30s: %s\n" % (k, v))
