import sys
import pandas as pd

"""
Rudimentary differ for clinvar_alleles.tsv.gz

Usage: python diff_clinvar_alleles.py \
    <clinvar_alleles.A.tsv.gz> \
    <clinvar_alleles.B.tsv.gz>
"""

print "comparing A: {} and B: {}".format(sys.argv[1], sys.argv[2])

sep = "#" * 25

print sep
INDEX = ['chrom', 'pos', 'ref', 'alt', 'allele_id']
df_a = pd.read_csv(sys.argv[1], sep="\t", compression='gzip',
                   index_col=False).sort_values(INDEX).set_index(INDEX)
print "A: {} rows".format(df_a.shape[0])
df_b = pd.read_csv(sys.argv[2], sep="\t", compression='gzip',
                   index_col=False).sort_values(INDEX).set_index(INDEX)
print "B: {} rows".format(df_b.shape[0])
rows_a = set(df_a.index)
rows_b = set(df_b.index)
common_rows = sorted(rows_a & rows_b)
print "{} rows in common".format(len(common_rows))
rows_a_not_b = list(rows_a - rows_b)
print "{} rows in A but not B, e.g.: {}".format(
    len(rows_a_not_b), ", ".join(map(str, rows_a_not_b[:5])))
rows_b_not_a = list(rows_b - rows_a)
print "{} rows in B but not A, e.g.: {}".format(
    len(rows_b_not_a), ", ".join(map(str, rows_b_not_a[:5])))
print sep

cols_a = set(df_a.columns)
cols_b = set(df_b.columns)
common_columns = sorted(cols_a & cols_b)
print "{} columns in common: {}".format(
    len(common_columns), ", ".join(common_columns))
print "{} columns in A but not B: {}".format(
    len(cols_a - cols_b), ", ".join(sorted(cols_a - cols_b)))
print "{} columns in B but not A: {}".format(
    len(cols_b - cols_a), ", ".join(sorted(cols_b - cols_a)))

print sep

df_a = df_a[common_columns].loc[common_rows]
df_b = df_b[common_columns].loc[common_rows]

for col in common_columns:
    print "comparing:", col
    matches = df_a[col] == df_b[col]
    matches |= df_a[col].isnull() & df_b[col].isnull()
    n_match = matches.sum()
    n_diff = matches.shape[0] - n_match
    print "{} rows match, {} are different".format(n_match, n_diff)
    if n_diff:
        print "example differences:"
        showed = 0
        for ix, is_match in matches.iteritems():
            if not is_match:
                print "\t{}: '{}' vs '{}'".format(
                    ix, df_a.loc[ix][col].values, df_b.loc[ix][col].values)
                showed += 1
                if showed >= 5:
                    break
    print sep