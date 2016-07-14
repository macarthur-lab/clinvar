"""
Script for generating a new clinvar table with the ExAC fields below added to each clinvar variant that's in ExAC
"""

import argparse
import pysam

NEEDED_EXAC_FIELDS = [
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
p.add_argument("-o", "--output-table", help="Output filename", default="clinvar_with_exac.tsv")
args = p.parse_args()

"""Clinvar table header:
chrom    pos    ref    alt    mut    measureset_id    symbol    clinical_significance    review_status    hgvs_c    hgvs_p    all_submitters    all_traits    all_pmids    pathogenic    conflicted
"""
exac_f = pysam.TabixFile(args.exac_sites_vcf)
clinvar_f = open(args.clinvar_table)
clinvar_header = next(clinvar_f).rstrip('\n').split('\t')
clinvar_with_exac_header = clinvar_header + NEEDED_EXAC_FIELDS
clinvar_with_exac_f = open(args.output_table, "w")
clinvar_with_exac_f.write("\t".join(clinvar_with_exac_header)+"\n")
for i, clinvar_row in enumerate(clinvar_f):
    clinvar_fields = clinvar_row.rstrip('\n').split('\t')
    clinvar_dict = dict(zip(clinvar_header, clinvar_fields))

    chrom = clinvar_dict['chrom']
    found_exac_values = False
    if chrom != 'MT':
        pos = int(clinvar_dict['pos'])

        # retrieve ExAC variant - pysam.fetch(..) sometimes returns more than 1 vcf record, so need to filter here
        matching_exac_rows = [r for r in exac_f.fetch(chrom, pos-1, pos) if r.startswith("%s\t%s\t" % (chrom, pos))]  

        assert len(matching_exac_rows) < 2

        if matching_exac_rows:
            exac_row_fields = matching_exac_rows[0].split('\t')
            filter_value = exac_row_fields[6]
            info_fields = [tuple(kv.split('=')) for kv in exac_row_fields[7].split(';')]
            info_fields = filter(lambda kv: kv[0] in NEEDED_EXAC_FIELDS_SET, info_fields)
            assert len(info_fields) == len(NEEDED_EXAC_FIELDS_SET)
            exac_column_values = [kv[1] for kv in info_fields]
            clinvar_with_exac_f.write("\t".join(clinvar_fields + exac_column_values) + "\n")
            found_exac_values = True

    if not found_exac_values:
        clinvar_with_exac_f.write("\t".join(clinvar_fields + EXAC_EMPTY_COLUMN_VALUES) + "\n")
