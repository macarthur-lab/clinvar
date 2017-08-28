#!/usr/bin/env python
import sys
import pandas as pd

from parse_clinvar_xml import HEADER

FINAL_HEADER = HEADER + ['gold_stars', 'conflicted']


def join_variant_summary_with_clinvar_alleles(
        variant_summary_table, clinvar_alleles_table,
        genome_build_id="GRCh37"):
    variant_summary = pd.read_csv(variant_summary_table, sep="\t",
                                  index_col=False, compression="gzip",low_memory=False)
    print "variant_summary raw", variant_summary.shape

    clinvar_alleles = pd.read_csv(clinvar_alleles_table, sep="\t",
                                  index_col=False, compression="gzip",low_memory=False)
    print "clinvar_alleles raw", clinvar_alleles.shape

    # use lowercase names and replace . with _ in column names:
    variant_summary = variant_summary.rename(columns=dict(
        (col, col.lower().replace(".", "_"))
        for col in variant_summary.columns))
    # rename first column to allele_id:
    variant_summary = variant_summary.rename(
        columns={variant_summary.columns[0]: "allele_id"})

    # extract relevant columns for the correct assembly and
    # rename clinicalsignificance, reviewstatus, lastevaluated:
    variant_summary = variant_summary[
        variant_summary['assembly'] == genome_build_id]
    variant_summary = variant_summary[
        ['allele_id', 'clinicalsignificance', 'reviewstatus','lastevaluated']]
    variant_summary = variant_summary.rename(
        columns={'clinicalsignificance': 'clinical_significance',
                 'reviewstatus': 'review_status',
                 'lastevaluated':'last_evaluated'})

    # remove the duplicated records in variant summary due to alternative loci such as PAR but would be problematic for rare cases like translocation
    variant_summary=variant_summary.drop_duplicates()
    print "variant_summary after filter", variant_summary.shape

    # remove clinical_significance and review_status from clinvar_alleles:
    clinvar_alleles = clinvar_alleles.drop(
        ['clinical_significance', 'review_status','last_evaluated'], axis=1)

    # pandas is sensitive to some rows having allele_id joined on ;, causing
    # an object dtype, with some entries being ints and others strs
    clinvar_alleles['allele_id'] = clinvar_alleles['allele_id'].astype(str)
    variant_summary['allele_id'] = variant_summary['allele_id'].astype(str)

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
        'practice guideline': 4,
        '-':'-'
    }
    df['gold_stars'] = df.review_status.map(gold_star_map)

    # The use of expressions on clinical significance on ClinVar aggregate records (RCV) https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#conflicts
    # conflicted = 1 if using "conflicting"
    df['conflicted'] = df['clinical_significance'].str.contains(
        "onflicting", case=False).astype(int)

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
