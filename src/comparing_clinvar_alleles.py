"""
Super simple script to diff two tables.
Lists out which columns are in common, which are not.
For each column, tallies up discordances and prints some examples.

Usage:
    python comparing_clinvar_alleles <table 1> <table 2> <index column>

The table can be a .tsv/.csv/.tsv.gz/.csv.gz

For clinvar_alleles, use variation_id
For clinvar_allele_trait_pairs, use rcv
"""
import sys
import pandas as pd

# quick & dirty -- could use a nicer CLI here
file_a = sys.argv[1]
file_b = sys.argv[2]
index_col = sys.argv[3]


def load_table(fname):
    sep = "\t" if ".tsv" in fname else ","
    compression = 'gzip' if fname.endswith('.gz') else none
    df = pd.read_csv(fname, index_col=False, compression=compression, sep=sep)
    n_unique = len(set(df[index_col]))
    assert n_unique == df.shape[0], \
        "Given index column {}, but there are {} unique values for {} rows".format(
            index_col, n_unique, df.shape[0])
    df = df.set_index(index_col)
    return df

df_a = load_table(file_a)
df_b = load_table(file_b)

cols_a = set(df_a.columns)
cols_b = set(df_b.columns)
cols_common = sorted(cols_a & cols_b)
cols_a_not_b = cols_a - cols_b
print("Columns in {} but not {}:".format(file_a, file_b))
print("\n".join(sorted(cols_a_not_b)))
print("")
cols_b_not_a = cols_b - cols_a
print("Columns in {} but not {}:".format(file_b, file_a))
print("\n".join(sorted(cols_b_not_a)))
print("")

rows_a = set(df_a.index)
rows_b = set(df_b.index)
rows_common = sorted(rows_a & rows_b)
n_common = len(rows_common)
print("Rows in common: {}".format(n_common))
rows_a_not_b = rows_a - rows_b
print("Rows in {} but not {}: {}".format(file_a, file_b, len(rows_a_not_b)))
rows_b_not_a = rows_b - rows_a
print("Rows in {} but not {}: {}".format(file_b, file_a, len(rows_b_not_a)))

df_a = df_a.loc[rows_common]
df_b = df_b.loc[rows_common]

for col in cols_common:
    print(col)
    exact_matches = df_a[col] == df_b[col]
    exact_matches |= pd.isnull(df_a[col]) & pd.isnull(df_b[col])  # force nan == nan
    n_exact = exact_matches.sum()
    print("\t{}/{} ~ {:,.3%} exact matches".format(n_exact, n_common, n_exact / float(n_common)))
    if n_exact < n_common:
        print("\texample differences:")
        for ix in df_a.index[~exact_matches][:5]:
            print("\t\t{}={}: '{}' vs '{}'".format(index_col, ix,
                                                   df_a.loc[ix, col], df_b.loc[ix, col]))
