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
    last_unique_id = ''
    last_data = {}
    counter = 0
    for line in input_fileobject.readlines():
        data = dict(zip(column_names,line.strip('\n').split('\t')))
        unique_id = '_'.join([data['chrom'], str(data['pos']), data['ref'], data['alt']])
        if counter == 0:
            pass
        elif unique_id == last_unique_id:
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
    print data1['worst_assertion']
    print data2['worst_assertion']
    if data1['mut'] == 'REF' and data2['mut'] == 'ALT':
        return data2
    elif data1['mut'] == 'ALT' and data2['mut'] == 'REF':
        return data1
    else:
        combined_data = data1 # this sets defaults, now we fix it:
    # concatenate catch-all fields
    for key in ['accession', 'all_traits', 'all_submitters', 'all_clnsig', 'all_pmids']:
        if data1[key] != data2[key]:
            combined_data[key] = ";".join([data1[key],data2[key]])
    # decide priority for unique fields
    # if one has higher review status, it wins...
    if revstat_ranking.get(data2['highest_revstat'],0) > revstat_ranking.get(data1['highest_revstat'],0):
        combined_data['highest_revstat'] = data2['highest_revstat']
        combined_data['worst_assertion'] = data2['worst_assertion']
    # if both have equal review status, we decide based on pathogenicity:
    elif revstat_ranking.get(data2['highest_revstat'],0) == revstat_ranking.get(data1['highest_revstat'],0):
        if clnsig_ranking_path_spectrum.get(data2['worst_assertion'],0) > clnsig_ranking_path_spectrum.get(data1['worst_assertion'],0):
            combined_data['worst_assertion'] = data2['worst_assertion']
    # if they are different or at least one is conflicted, then it is conflicted
    if data2['worst_assertion'] != data1['worst_assertion'] or data2['conflicted'] == '1' or data1['conflicted'] == '1':
        combined_data['conflicted'] = '1'
    elif data2['highest_revstat'] == 'classified by single submitter' and data1['highest_revstat'] == 'classified by single submitter':
        # only if there are two single submitter classifications and they are NOT conflicted, then call it "multiple submitters"
        combined_data['highest_revstat'] = 'classified by multiple submitters'
    return combined_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='De-duplicate the output from parse_clinvar_xml.py')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
                       default=sys.stdin)
    args = parser.parse_args()
    dedup_clinvar(args.infile)


