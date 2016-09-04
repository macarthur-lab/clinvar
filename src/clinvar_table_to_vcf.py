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

    print("""
##fileformat=VCFv4.1
##source=clinvar
##INFO=<ID=MUT,Number=1,Type=String,Description="MUT">
##INFO=<ID=MEASURESET_ID,Number=1,Type=String,Description="MEASURESET_ID">
##INFO=<ID=SYMBOL,Number=1,Type=String,Description="SYMBOL">
##INFO=<ID=CLINICAL_SIGNIFICANCE,Number=1,Type=String,Description="CLINICAL_SIGNIFICANCE">
##INFO=<ID=REVIEW_STATUS,Number=1,Type=String,Description="REVIEW_STATUS">
##INFO=<ID=ALL_SUBMITTERS,Number=1,Type=String,Description="ALL_SUBMITTERS">
##INFO=<ID=ALL_TRAITS,Number=1,Type=String,Description="ALL_TRAITS">
##INFO=<ID=ALL_PMIDS,Number=1,Type=String,Description="ALL_PMIDS">
##INFO=<ID=INHERITANCE_MODES,Number=1,Type=String,Description="INHERITANCE_MODES">
##INFO=<ID=AGE_OF_ONSET,Number=1,Type=String,Description="AGE_OF_ONSET">
##INFO=<ID=PREVALENCE,Number=1,Type=String,Description="PREVALENCE">
##INFO=<ID=DISEASE_MECHANISM,Number=1,Type=String,Description="DISEASE_MECHANISM">
##INFO=<ID=ALL_PATHOGENIC,Number=1,Type=String,Description="ALL_PATHOGENIC">
##INFO=<ID=ALL_CONFLICTED,Number=1,Type=String,Description="ALL_CONFLICTED">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=MT,length=16569>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##reference=Homo_sapiens_assembly19.fasta
""".strip())

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
                    "inheritance_modes", "age_of_onset", "prevalence", 'disease_mechanism',
                    "pathogenic", "conflicted"]:
            if pd.isnull(table_row[key]):
                continue
            value = str(table_row[key])
            value = re.sub('\s*[,]\s*', '..', value)  # replace , with ..
            value = re.sub('\s*[;]\s*', '|', value)  # replace ; with |
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
