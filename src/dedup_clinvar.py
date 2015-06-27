#!/usr/bin/env python

import re
import sys
import gzip
import argparse

from parse_clinvar_xml import *

# recommended usage:
# ./dedup_clinvar.py < clinvar_table_sorted.tsv > clinvar_table_dedup.tsv

'''
Run through a sorted clinvar_table.tsv file from the parse_clinvar_xml script, and make it unique on CHROM POS REF ALT
'''
def dedup_clinvar(input_fileobject):
    header = input_fileobject.readline()
    sys.stdout.write(header)
    column_names = header.strip('\n').split('\t')
    first_data_line = input_fileobject.readline()
    last_data = dict(zip(column_names,first_data_line.strip('\n').split('\t')))
    last_unique_id = '_'.join([last_data['chrom'], str(last_data['pos']), last_data['ref'], last_data['alt']])
    counter = 0
    for line in input_fileobject.readlines():
        data = dict(zip(column_names,line.strip('\n').split('\t')))
        unique_id = '_'.join([data['chrom'], str(data['pos']), data['ref'], data['alt']])
        if unique_id == last_unique_id:
            data = dedup_records(data, last_data)
        else:
            # note that using a comprehension instead of just data.values() preserves column order
            sys.stdout.write('\t'.join([last_data[colname] for colname in column_names])+'\n')
        last_data = data
        last_unique_id = unique_id
        counter += 1
    sys.stdout.write('\t'.join([last_data[colname] for colname in column_names])+'\n')

'''
De-duplicate two ClinVar records
'''
def dedup_records(data1, data2):
    # if one is REF and one is ALT, discard the REF
    if data1['mut'] == 'REF' and data2['mut'] == 'ALT':
        return data2
    elif data1['mut'] == 'ALT' and data2['mut'] == 'REF':
        return data1
    else:
        combined_data = data1 # this sets defaults, now we fix it:
    # discard one MeasureSet ID if there are two
    if data1['measureset_id'] != data2['measureset_id']:
        combined_data['measureset_id'] = str(min(map(int,[data1['measureset_id'],data2['measureset_id']])))
    combined_data['all_pmids'] = ','.join(set(data1['all_pmids'].split(',') + data2['all_pmids'].split(',')))
    combined_data['all_submitters'] = ';'.join(set(data1['all_submitters'].split(';') + data2['all_submitters'].split(';')))
    return combined_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='De-duplicate the output from parse_clinvar_xml.py')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
                       default=sys.stdin)
    args = parser.parse_args()
    dedup_clinvar(args.infile)


