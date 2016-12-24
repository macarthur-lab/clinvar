#!/usr/bin/env python

import re
import sys
import gzip
import argparse

from parse_clinvar_xml import *

# recommended usage:
# ./dedup_clinvar.py < clinvar_table_sorted.tsv > clinvar_table_dedup.tsv

#modified by Xioalei Zhang
'''
Run through a sorted clinvar_table.tsv file from the parse_clinvar_xml script, and make it unique on CHROM POS REF ALT AlleleID
'''
def dedup_clinvar(input_fileobject,outfile):
    header = input_fileobject.readline()
    outfile.write(header)
    column_names = header.strip('\n').split('\t')
    first_data_line = input_fileobject.readline()
    last_data = dict(zip(column_names,first_data_line.strip('\n').split('\t')))
    last_unique_id = '_'.join([last_data['chrom'], str(last_data['pos']), last_data['ref'], last_data['alt'],last_data['allele_id']])    
    counter = 0
    for line in input_fileobject.readlines():
        data = dict(zip(column_names,line.strip('\n').split('\t')))
        unique_id = '_'.join([data['chrom'], str(data['pos']), data['ref'], data['alt'],data['allele_id']])
        if unique_id == last_unique_id:
            data = dedup_records(data, last_data)
        else:
            # note that using a comprehension instead of just data.values() preserves column order
            # the next line (data) is not duplicated as the current line(last_data) then just print last_data
            outfile.write('\t'.join([last_data[colname] for colname in column_names])+'\n')
        last_data = data
        last_unique_id = unique_id
        counter += 1
    outfile.write('\t'.join([last_data[colname] for colname in column_names])+'\n')

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
    
    # concatenate columns that may have lists of values    
    for column_name in ('measureset_type','measureset_id','rcv',
        'symbol','hgvs_c','hgvs_p','molecular_consequence','clinical_significance', 
        'pathogenic', 'benign', 'conflicted', 'review_status', 'gold_stars','all_submitters',
        'all_traits','all_pmids', 'inheritance_modes', 'age_of_onset','prevalence', 'disease_mechanism', 
        'origin', 'xrefs'):
        combined_data[column_name] = ';'.join(set(filter(lambda s: s,
                data1[column_name].split(';') + 
                data2[column_name].split(';'))))

    return combined_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='De-duplicate the output from parse_clinvar_xml.py')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
                       default=sys.stdin)
    parser.add_argument('-o','--outfile',nargs='?',type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()
    dedup_clinvar(args.infile,args.outfile)
