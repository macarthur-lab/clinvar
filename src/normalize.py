#!/usr/bin/env python

'''
This script is a python implementation of the algorithm for variant
normalization described by Tan et al 2015:

Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants.
Bioinformatics. 2015 Jul 1;31(13):2202-4. doi: 10.1093/bioinformatics/btv112.
Epub 2015 Feb 19. PubMed PMID: 25701572.

The authors have made a C++ implementation available in vt as vt normalize
And the source code is viewable here: https://github.com/atks/vt

For our purposes, we wanted a Python implementation so that we could
build end-to-end scripts in Python.

If you use this, please cite Tan et al 2015.

A note about when this is useful. In VCFs generated with GATK (or probably
other tools) from short read sequencing data, variants are already left-aligned
but may be non-minimal to the extent that indels overlap with other variants.
For those cases, minimal_representation.py is sufficient to convert variants
to minimal representation. However, genomic coordinates converted from HGVS
(we have encountered this when parsing the ClinVar XML dump) may be not only
non-minimal but also right-aligned rather than left-aligned, and may contain
hyphens. For those situations, use this script (or just run vt normalize).

Usage: normalize.py -R $b37ref < bad_file.txt > good_file.txt
'''

import sys
import pysam
import argparse

'''
An Error class for rare cases where REF == ALT (seen in ClinVar XML)
'''
class RefEqualsAltError(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)

'''
An Error class for REF or ALT values that are not valid nucleotide sequences
'''
class InvalidNucleotideSequenceError(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)

'''
An Error class for variants where the REF does not match the reference genome
'''
class WrongRefError(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)
   

'''
Accepts a pysam FastaFile object pointing to the reference genome, and
chrom, pos, ref, alt genomic coordinates, and normalizes them.
'''
def normalize(pysam_fasta, chrom, pos, ref, alt):
    pos = int(pos) # make sure position is an integer
    ref = ref.upper()
    alt = alt.upper()
    # Remove variants that contain invalid nucleotides
    if any(letter not in ['A','C','G','T','N','-'] for letter in ref + alt):
        raise InvalidNucleotideSequenceError('Invalid nucleotide sequence: %s %s %s %s'%(chrom, pos, ref, alt))
    # use blanks instead of hyphens
    if ref == '-':
        ref = ''
    if alt == '-':
        alt = ''
    # check whether the REF is correct
    true_ref = pysam_fasta.fetch(chrom, pos - 1, pos - 1 + len(ref))
    if ref != true_ref:
        raise WrongRefError('Incorrect REF value: %s %s %s %s (actual REF should be %s)'%(chrom, pos, ref, alt, true_ref))
    # Prevent infinte loops in cases where REF == ALT.
    # We have encountered this error in genomic coordinates from the ClinVar XML file
    if ref == alt:
        raise RefEqualsAltError('The REF and ALT allele are the same: %s %s %s %s'%(chrom, pos, ref, alt))
    # Time-saving shortcut for SNPs that are already minimally represented
    if len(ref) == 1 and len(alt) == 1 and ref in ['A','C','G','T'] and alt in ['A','C','G','T']:
        return chrom, pos, ref, alt
    # This first loop left-aligns and removes excess nucleotides on the right.
    # This is Algorithm 1 lines 1-6 from Tan et al 2015
    keep_working = True
    while keep_working:
        keep_working = False
        if len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            keep_working = True
        if len(ref) == 0 or len(alt) == 0:
            preceding_base = pysam_fasta.fetch(chrom, pos-2, pos-1)
            ref = preceding_base + ref
            alt = preceding_base + alt
            pos = pos - 1
            keep_working = True
    # This second loop removes excess nucleotides on the left. This is Algorithm 1 lines 7-8.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    return chrom, pos, ref, alt

'''
This function takes a tab-delimited file with a header line containing columns
named chrom, pos, ref, and alt, plus any other columns. It normalizes the
chrom, pos, ref, and alt, and writes all columns out to another file.
'''
def normalize_tab_delimited_file(infile, outfile, reference_fasta, verbose=True):
    pysam_fasta = pysam.FastaFile(reference_fasta) # create a pysam object of the reference genome
    header = infile.readline() # get header of input file
    columns = header.strip('\n').split('\t')  # parse col names 
    outfile.write('\t'.join(columns) + '\n') # write header line plus the CpG col to be generated
    counter = 0
    ref_equals_alt = 0
    wrong_ref = 0
    invalid_nucleotide = 0
    for line in infile.readlines():
        data = dict(zip(columns,line.strip('\n').split('\t')))
        # fill the data with blanks for any missing data
        for column in columns:
            if column not in data.keys():
                data[column] = ''
        pos = int(data['pos'])
        try:
            data['chrom'], pos, data['ref'], data['alt'] = normalize(pysam_fasta, data['chrom'], pos, data['ref'], data['alt'])
        except RefEqualsAltError as e:
            sys.stderr.write('\n'+str(e)+'\n')
            ref_equals_alt += 1
            continue
        except WrongRefError as e:
            sys.stderr.write('\n'+str(e)+'\n')
            wrong_ref += 1
            continue
        except InvalidNucleotideSequenceError as e:
            sys.stderr.write('\n'+str(e)+'\n')
            invalid_nucleotide += 1
            continue
        data['pos'] = str(pos)
        outfile.write('\t'.join([data[column] for column in columns]) + '\n')
        counter += 1
        if verbose:
            sys.stderr.write("\r%s records processed\n"%(counter))
    outfile.write('\n\n')
    if verbose:
        sys.stderr.write("Final counts of variants discarded:\nREF == ALT: %s\nWrong REF: %s\nInvalid nucleotide: %s\n"%(ref_equals_alt, wrong_ref, invalid_nucleotide))

'''
Battery of test cases for normalize
'''
def test_normalize(pysam_fasta):
    sys.stdout.write(str(normalize(pysam_fasta, '7', 117199646, 'CTT', '-'))+'\n') # HGVS translation of CFTR p.F508del, should be ('7', 117199644, 'ATCT', 'A')
    sys.stdout.write(str(normalize(pysam_fasta, '13', 32914438, 'T', '-'))+'\n') # HGVS translation of a BRCA2 Ashkenazi founder variant, should be ('13', 32914437, 'GT', 'G')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python implementation of vt normalize')
    parser.add_argument('-R', '--reference_fasta', type=str, default='',
        help="Path to FASTA file of reference genome. Must be samtools faidx'ed")
    args = parser.parse_args()
    normalize_tab_delimited_file(sys.stdin, sys.stdout, args.reference_fasta)
