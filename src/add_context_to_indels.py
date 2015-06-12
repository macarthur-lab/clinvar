#!/usr/bin/env python

# Genomic coordinates that have been parsed from HGVS notation often lack the one base of 
# reference genome context that is required by the VCF spec.
# For instance, the ClinVar XML release refers to an Ashkenazi Jewish BRCA2 founder variant as:
# 13    32914438    T   -
# When the VCF spec-compliant representation would be:
# 13    32914438    GT   G
# Fixing these problems requires looking up the reference base (in this case G), which we do using pysam

# Usage: add_context_to_indels.py -R $b37ref < bad_file.txt > good_file.txt

import sys
import pysam
import argparse

def add_context(infile,outfile,reference_fasta):
    ref = pysam.FastaFile(reference_fasta) # create a pysam object of the reference genome
    header = infile.readline() # get header of input file
    columns = header.strip('\n').split('\t')  # parse col names 
    outfile.write('\t'.join(columns) + '\n') # write header line plus the CpG col to be generated
    for line in infile.readlines():
        data = dict(zip(columns,line.strip('\n').split('\t')))
        # fill the data with blanks for any missing data
        for column in columns:
            if column not in data.keys():
                data[column] = ''
        preceding_base = ref.fetch(data['chrom'], int(data['pos'])-2, int(data['pos'])-1)
        if data['ref'] == '-' or data['alt'] == '-': # if it's an indel with missing context
            data['pos'] = str(int(data['pos'])-1) # reduce position by 1 base
            data['ref'] = preceding_base + data['ref'].replace('-','') # prepend preceding base, remove '-'
            data['alt'] = preceding_base + data['alt'].replace('-','') # prepend preceding base, remove '-'
        outfile.write('\t'.join([data[column] for column in columns]) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add indel context to fix genomic coordinates that violate the VCF spec')
    parser.add_argument('-R', '--reference_fasta', type=str, default='',
        help="Path to FASTA file of reference genome. Must be samtools faidx'ed")
    args = parser.parse_args()
    add_context(sys.stdin, sys.stdout, args.reference_fasta)
