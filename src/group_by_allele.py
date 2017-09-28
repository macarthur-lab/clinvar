#!/usr/bin/env python

import argparse
import gzip
import sys

from parse_clinvar_xml import HEADER
# recommended usage:
# ./group_by_allele.py < clinvar_combined.tsv > clinvar_alleles.tsv


def group_by_allele(infile, outfile):
    """Run through a sorted clinvar_table.tsv file from the parse_clinvar_xml script, and make it unique on CHROM POS REF ALT

    Args:
        infile: Input file stream for reading clinvar_table_sorted.tsv
        outfile: Output file stream to write to.
    """

    header = next(infile)
    outfile.write(header)
    column_names = header.strip('\n').split('\t')

    last_data = None
    last_unique_id = None
    counter = 0

    for line in infile:
        data = dict(zip(column_names, line.strip('\n').split('\t')))
        unique_id = '_'.join([data['chrom'], str(data['pos']), data['ref'], data['alt']])
        if unique_id == last_unique_id:
            data = group_alleles(last_data, data)
        elif last_data is not None:
            # note that using a comprehension instead of just data.values() preserves column order
            # the next line (data) is not duplicated as the current line(last_data) then just print last_data
            outfile.write('\t'.join([last_data[colname] for colname in column_names])+'\n')
        last_data = data
        last_unique_id = unique_id
        counter += 1

    if last_data is not None:
        outfile.write('\t'.join([last_data[colname] for colname in column_names])+'\n')
    else:
        raise ValueError("%s has 0 records" % infile)

def group_alleles(data1, data2):
    """Group two variants with same genomic coordinates.

    Args:
        data1: dictionary of column-name, value pairs for table record #1
        data2: dictionary of column-name, value pairs for table record #2
    """

    if (data1['chrom'], data1['pos'], data1['ref'], data1['alt']) != (data2['chrom'], data2['pos'], data2['ref'], data2['alt']):
        raise ValueError("data1 variant id != data2 variant id: %s != %s" % (data1, data2))

    combined_data = data1  # this sets defaults, now we fix it:

    # 'pathogenic', 'benign', 'conflicted', 'gold_stars',
    # concatenate columns that may have lists of values
    loc_column = ['chrom','pos','ref','alt']
    num_field = ['pathogenic', 'likely_pathogenic','uncertain_significance','likely_benign', 'benign']
    info_column = [x for x in HEADER if x not in loc_column and x not in num_field]

    for column_name in info_column:
        all_non_empty_values = filter(lambda s: s, data1[column_name].split(';') + data2[column_name].split(';'))
        # deduplicate values, while preserving order
        deduplicated_values = []
        for value in all_non_empty_values:
            if value not in deduplicated_values:
                deduplicated_values.append(value)

        combined_data[column_name] = ';'.join(deduplicated_values)

    for column_name in num_field:
        combined_data[column_name]=str(int(data1[column_name])+int(data2[column_name]))

    return combined_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='De-duplicate the output from parse_clinvar_xml.py')
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    if args.infile.name.endswith(".gz"):
        args.infile.close()
        args.infile = gzip.open(args.infile.name)
    
    if args.outfile.name.endswith(".gz"):
        args.outfile.close()
        args.outfile = gzip.open(args.outfile.name, 'w')
    
    group_by_allele(args.infile, args.outfile)
