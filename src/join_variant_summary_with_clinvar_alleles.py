#!/usr/bin/env python
import sys
import pandas as pd

from parse_clinvar_xml import HEADER

FINAL_HEADER = HEADER + ['gold_stars', 'pathogenic', 'benign', 'conflicted',
                         'uncertain']


def join_variant_summary_with_clinvar_alleles(
        variant_summary_table, clinvar_alleles_table,
        genome_build_id="GRCh37"):
    variant_summary = pd.read_csv(variant_summary_table, sep="\t",
                                  index_col=False, compression="gzip")
    print "variant_summary raw", variant_summary.shape

    clinvar_alleles = pd.read_csv(clinvar_alleles_table, sep="\t",
                                  index_col=False, compression="gzip")
    print "clinvar_alleles raw", clinvar_alleles.shape

    # use lowercase names and replace . with _ in column names:
    variant_summary = variant_summary.rename(columns=dict(
        (col, col.lower().replace(".", "_"))
        for col in variant_summary.columns))
    # rename first column to allele_id:
    variant_summary = variant_summary.rename(
        columns={variant_summary.columns[0]: "allele_id"})

    # extract relevant columns for the correct assembly and
    # rename clinicalsignificance, reviewstatus:
    variant_summary = variant_summary[
        variant_summary['assembly'] == genome_build_id]
    variant_summary = variant_summary[
        ['allele_id', 'clinicalsignificance', 'reviewstatus']]
    variant_summary = variant_summary.rename(
        columns={'clinicalsignificance': 'clinical_significance',
                 'reviewstatus': 'review_status'})
    print "variant_summary after filter", variant_summary.shape

    # remove clinical_significance and review_status from clinvar_alleles:
    clinvar_alleles = clinvar_alleles.drop(
        ['clinical_significance', 'review_status'], axis=1)

    # remove rows with multiple allele_ids
    clinvar_alleles = clinvar_alleles[
        ~(clinvar_alleles['allele_id'].str.contains(";") == True)]
    # have to do it this way because allele_id gets a dtype of object, being
    # a mix of ints and strs; alternate way would be to cast
    # variant_summary.allele_id to str and let pd.merge handle it
    clinvar_alleles['allele_id'] = clinvar_alleles.allele_id.astype(int)

    print "clinvar_alleles after filter", clinvar_alleles.shape

    # join the two tables on allele_id:
    df = pd.merge(clinvar_alleles, variant_summary,
                  on='allele_id', how='inner')
    print "merged raw", df.shape

    # map review_status to gold starts:
    gold_star_map = {
        'no assertion provided': 0,
        'no assertion for the individual variant': 0,
        'no assertion criteria provided': 0,
        'criteria provided, single submitter': 1,
        'criteria provided, conflicting interpretations': 1,
        'criteria provided, multiple submitters, no conflicts': 2,
        'reviewed by expert panel': 3,
        'practice guideline': 4
    }
    df['gold_stars'] = df.review_status.map(gold_star_map)

    # pathogenic = 1 if at least one submission says pathogenic or likely
    # pathogenic, 0 otherwise
    df['pathogenic'] = df['clinical_significance'].str.contains(
        "pathogenic", case=False)
    # benign = 1 if at least one submission says benign or likely benign,
    # 0 otherwise
    df['benign'] = df['clinical_significance'].str.contains(
        "benign", case=False)
    # conflicted = 1 if pathogenic == 1 and benign == 1 or if the significance
    # string contains "conflicting data"
    df['conflicted'] = df['pathogenic'] & df['benign']
    df['conflicted'] |= df['clinical_significance'].str.contains(
        "conflicting data", case=False)
    # uncertain = 1 if the variant is of uncertain significance or if it is
    # conflicted
    df['uncertain'] = df['clinical_significance'].str.contains(
        "uncertain", case=False) | df['conflicted']
    for flag in ('pathogenic', 'benign', 'conflicted', 'uncertain'):
        df[flag] = df[flag].astype(int)

    # reorder columns just in case
    df = df.ix[:, FINAL_HEADER]
    print "merged final", df.shape

    return df


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python join_variant_summary_with_clinvar_alleles.py "
              "<variant_summary.txt.gz> "
              "<clinvar_alleles_grouped.tsv.gz> "
              "<out_name: clinvar_alleles_combined.tsv.gz> "
              "<genome build, e.g. GRCh37>")
        exit()
    variant_summary_table = sys.argv[1]
    clinvar_alleles_table = sys.argv[2]
    out_name = sys.argv[3]
    genome_build_id = sys.argv[4]
    assert out_name.endswith('.gz'), ("Provide a filename with .gz extension "
                                      "as the output will be gzipped")
    df = join_variant_summary_with_clinvar_alleles(
        variant_summary_table, clinvar_alleles_table, genome_build_id)
    df.to_csv(out_name, sep="\t", index=False, compression='gzip')
