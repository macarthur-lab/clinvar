"""
Basic consistency checks on the final clinvar table.
"""

import argparse
import gzip
import os

from pprint import pprint

p = argparse.ArgumentParser(description="Basic consistency checks on the final clinvar table")
p.add_argument("clinvar_table")
args = p.parse_args()

clinvar_table = args.clinvar_table
if not os.path.isfile(clinvar_table):
    p.error("%s doesn't exist" % clinvar_table)

CHROMS = list(map(str, range(1, 23))) + ['X', 'Y', 'MT']

f = gzip.open(clinvar_table) if clinvar_table.endswith('gz') else open(clinvar_table)
header = next(f).strip('\n').split('\t')
#print(header)
#print(next(f).strip('\n').split('\t'))
errors_counter = 0
for i, line in enumerate(f):
    #print(line.strip('\n').split('\t'))

    record = dict(zip(header, line.strip('\n').split('\t')))
    try:
        assert record['chrom'] in CHROMS, 'Unexpected "chrom" column value: ' + record['chrom']
        assert int(record['pos']) > 0 and int(record['pos']) < 3*10**8, 'Unexpected "pos" column value: ' + record['pos']
        assert all(b in 'ACGTN' for b in record['ref']), 'Unexpected "ref" column value: ' + record['ref']
        assert all(b in 'ACGTN' for b in record['alt']), 'Unexpected "alt" column value: ' + record['alt']  # there's one clinvar record with ALT = "NTGT". Not sure how to handle it.
        assert record['measureset_type'] in ["Variant", "Haplotype", "CompoundHeterozygote", "Phase unknown", "Distinct chromosomes", "CompoundHeterozygote;Haplotype", " Variant;gene-variant"], \
            'Unexpected "measureset_type" column value: ' + record['measureset_type']  # there's one clinvar record with ALT = "NTGT". Not sure how to handle it.

        assert len(map(int, record['measureset_id'].split(';'))) > 0, 'Unexpected "measureset_id" column value: ' + record['measureset_id']
        assert len(map(lambda rcv: int(rcv.strip('RCV')), record['rcv'].split(';'))) > 0, 'Unexpected "rcv" column value: ' + record['rcv']
        assert int(record['allele_id']) > 0, 'Unexpected "rcv" column value: ' + record['allele_id']
        assert len(record['hgvs_c']) == 0 or "c." in record['hgvs_c'], 'Unexpected "hgvs_c" column value: ' + record['hgvs_c']
        assert len(record['hgvs_p']) == 0 or "p." in record['hgvs_p'], 'Unexpected "hgvs_p" column value: ' + record['hgvs_p']
        #assert record['molecular_consequence'], 'Unexpected "molecular_consequence" column value: ' + record['molecular_consequence']


        assert record['benign'] in ('0', '1'), 'Unexpected "benign" column value: ' + record['benign']
        assert record['pathogenic'] in ('0', '1'), 'Unexpected "pathogenic" column value: ' + record['pathogenic']
        assert record['conflicted'] in ('0', '1'), 'Unexpected "conflicted" column value: ' + record['conflicted']

        assert record['pathogenic'] == '1' and record['benign'] == '1' if record['conflicted'] == '1' else 1, 'Unexpected "conflicted" column value: ' + record['conflicted']

    except AssertionError as e:
        print("====================================")
        print("ERROR in %s - line %s: " % (clinvar_table, i))
        print(e)
        pprint(record)
        errors_counter += 1
if errors_counter > 0:
    p.error("%s errors found" % errors_counter)
