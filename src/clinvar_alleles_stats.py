#!/usr/bin/env python
import pandas as pd
import sys

"""
Summarizes some of the columns of clinvar_alleles.tsv.gz file
Usage: python clinvar_alleles_stats.py <clinvar_alleles.tsv.gz>
"""

alleles_name = sys.argv[1]

columns_to_summarize = [
    'variation_type', 'clinical_significance',
    'review_status', 'gold_stars', 'all_submitters',
    'inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism',
    'origin']

df = pd.read_csv(alleles_name, sep="\t", compression='gzip', index_col=False)

sep = "=" * 16

print "Columns:",
for i, col in enumerate(df.columns, start=1):
    print "{}: {},".format(i, col),
print

print sep

print "Total rows:", df.shape[0]

all_columns = list(df.columns)

for col in columns_to_summarize:
    print sep
    i = all_columns.index(col)
    print "column {}: {}".format(i + 1, col)
    print df[col].value_counts()

