"""
Basic consistency checks on the final clinvar table.
"""

import argparse
import gzip
import os

from pprint import pprint

p = argparse.ArgumentParser(description="Basic consistency checks on the final clinvar table")
p.add_argument("alleles_table_path")
args = p.parse_args()

alleles_table_path = args.alleles_table_path
if not os.path.isfile(alleles_table_path):
    p.error("%s doesn't exist" % alleles_table_path)

CHROMS = list(map(str, range(1, 23))) + ['X', 'Y', 'MT']

f = gzip.open(alleles_table_path) if alleles_table_path.endswith('gz') else open(alleles_table_path)
header = next(f).strip('\n').split('\t')
#print(header)
#print(next(f).strip('\n').split('\t'))
counter = 0
errors_counter = 0
for i, line in enumerate(f):
    counter += 1
    #print(line.strip('\n').split('\t'))

    record = dict(zip(header, line.strip('\n').split('\t')))
    try:
        assert record['chrom'] in CHROMS, 'Unexpected "chrom" column value: ' + record['chrom']
        assert int(record['pos']) > 0 and int(record['pos']) < 3*10**8, 'Unexpected "pos" column value: ' + record['pos']
        assert all(b in 'ACGTN' for b in record['ref']), 'Unexpected "ref" column value: ' + record['ref']
        assert all(b in 'ACGTN' for b in record['alt']), 'Unexpected "alt" column value: ' + record['alt']  # there's one clinvar record with ALT = "NTGT". Not sure how to handle it.
        assert record['variation_type'] in ["Variant", "Haplotype", "CompoundHeterozygote", "Phase unknown", "Distinct chromosomes", "CompoundHeterozygote;Haplotype", "Variant;gene-variant"], \
            'Unexpected "variation_type" column value: ' + record['variation_type']  # there's one clinvar record with ALT = "NTGT". Not sure how to handle it.

        assert len(map(int, record['variation_id'].split(';'))) > 0, 'Unexpected "variation_id" column value: ' + record['variation_id']
        assert len(map(lambda rcv: int(rcv.strip('RCV')), record['rcv'].split(';'))) > 0, 'Unexpected "rcv" column value: ' + record['rcv']
        assert int(record['allele_id']) > 0, 'Unexpected "rcv" column value: ' + record['allele_id']
        assert len(record['hgvs_c']) == 0 or "c." in record['hgvs_c'], 'Unexpected "hgvs_c" column value: ' + record['hgvs_c']
        assert len(record['hgvs_p']) == 0 or "p." in record['hgvs_p'], 'Unexpected "hgvs_p" column value: ' + record['hgvs_p']
        #assert record['molecular_consequence'], 'Unexpected "molecular_consequence" column value: ' + record['molecular_consequence']


    except AssertionError as e:
        print("====================================")
        print("ERROR in %s - line %s: " % (alleles_table_path, i))
        print(e)
        pprint(record)
        errors_counter += 1

assert ("multi" in alleles_table_path and counter > 100) or ("single" in alleles_table_path and counter > 10000), 'Table %s has only %s records' % (alleles_table_path, counter)
    
    

if errors_counter > 0:
    p.error("%s errors found" % errors_counter)
