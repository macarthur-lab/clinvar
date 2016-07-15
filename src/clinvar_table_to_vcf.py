import argparse
import collections
import gzip
import os
import re
import pandas as pd
import sys

def gzopen(path, mode, verbose=True):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)


def table_to_vcf(input_table_path):
    # validate args
    if not os.path.isfile(input_table_path):
        sys.exit("ERROR: %s not found" % input_table_path)

    # read input table. low_memory allows dtypes to be inferred
    t = pd.read_table(input_table_path, low_memory=False)

    missing_columns = {"chrom", "pos", "ref", "alt"} - set(t.columns)
    if missing_columns:
        sys.exit("ERROR: %s is missing columns: %s" % (input_table_path, str(missing_columns)))

    print("""##source=ClinVar
##INFO=<ID=MUT,Number=1,Type=String,Description="MUT">
##INFO=<ID=MUT,Number=1,Type=String,Description="MEASURESET_ID">
##INFO=<ID=MUT,Number=1,Type=String,Description="SYMBOL">
##INFO=<ID=MUT,Number=1,Type=String,Description="CLINICAL_SIGNIFICANCE">
##INFO=<ID=MUT,Number=1,Type=String,Description="REVIEW_STATUS">
##INFO=<ID=MUT,Number=1,Type=String,Description="ALL_SUBMITTERS">
##INFO=<ID=MUT,Number=1,Type=String,Description="ALL_TRAITS">
##INFO=<ID=MUT,Number=1,Type=String,Description="ALL_PMIDS">
##INFO=<ID=MUT,Number=1,Type=String,Description="ALL_PATHOGENIC">
##INFO=<ID=MUT,Number=1,Type=String,Description="ALL_CONFLICTED">""")
    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))
    for i, table_row in t.iterrows():
        vcf_row = []
        vcf_row.append(table_row["chrom"])
        vcf_row.append(table_row["pos"])
        vcf_row.append('.')  # ID
        vcf_row.append(table_row["ref"])
        vcf_row.append(table_row["alt"])
        vcf_row.append('.')  # QUAL
        vcf_row.append('.')  # FILTER

        info_field = collections.OrderedDict()

        # from VCF spec:
        #    INFO - additional information: (String, no white-space, semi-colons, or equals-signs permitted; commas are
        #    permitted only as delimiters for lists of values) INFO fields are encoded as a semicolon-separated series of short
        #    keys with optional values in the format: <key>=<data>[,data].
        for key in ["mut", "measureset_id", "symbol",
                    "clinical_significance", "review_status",
                    "all_submitters", "all_traits", "all_pmids",
                    "pathogenic", "conflicted"]:
            if pd.isnull(table_row[key]):
                continue
            value = str(table_row[key])
            value = re.sub('\s*[,;]\s*', '|', value)  # replace , or ; with |
            value = value.replace("=", " eq ").replace(" ", "_")
            
            info_field[key.upper()] = value
        vcf_row.append(";".join([key+"="+value for key, value in info_field.items()]))

        print("\t".join(map(str, vcf_row)))

    sys.stderr.write("Done\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_table_path', help="Tab-delimited input table")
    args = parser.parse_args()

    table_to_vcf(args.input_table_path)
