#!/usr/bin/env python

import re
import sys
import gzip
import argparse
from collections import defaultdict
import xml.etree.ElementTree as ET

# to test: ./parse_clinvar_xml.py -x clinvar_test.xml
# to run for reals: bsub -q priority -R rusage[mem=32] -oo cvxml.o -eo cvxml.e -J cvxml "./parse_clinvar_xml.py -x ClinVarFullRelease_00-latest.xml.gz -o clinvar_table.tsv"
# then sort it: cat clinvar_table.tsv | head -1 > clinvar_table_sorted.tsv; cat clinvar_table.tsv | tail -n +2 | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv

mentions_pubmed_regex = '(?:PubMed|PMID)(.*)' # group(1) will be all the text after the word PubMed or PMID
extract_pubmed_id_regex = '[^0-9]+([0-9]+)[^0-9](.*)' # group(1) will be the first PubMed ID, group(2) will be all remaining text

def parse_clinvar_tree(handle,dest=sys.stdout,verbose=True,mode='collapsed'):
    # print a header row
    dest.write(('\t'.join( ['chrom', 'pos', 'ref', 'alt', 'mut', 'measureset_id', 'all_submitters', 'all_traits', 'all_pmids'] ) + '\n').encode('utf-8'))
    counter = 0
    skipped_counter = defaultdict(int)
    for event, elem in ET.iterparse(handle):
        if event == 'end' and elem.tag == 'ClinVarSet':
            # find the GRCh37 VCF representation
            sequence_locations = elem.findall('.//SequenceLocation')
            grch37 = None
            for sequence_location in sequence_locations:
                if sequence_location.attrib.get('Assembly') == 'GRCh37':
                    if all(entry is not None for entry in [sequence_location.attrib.get(key) for key in ['referenceAllele','alternateAllele','start','Chr']]):
                        grch37 = sequence_location
            if grch37 is None:
                skipped_counter['missing SequenceLocation'] += 1
                elem.clear()
                continue # don't bother with variants that don't have a VCF location
            else:
                chrom = grch37.attrib['Chr']
                pos = grch37.attrib['start']
                ref = grch37.attrib['referenceAllele']
                alt = grch37.attrib['alternateAllele']
            measureset = elem.findall('.//MeasureSet')
            if measureset is None:
                skipped_counter['missing MeasureSet'] += 1
                elem.clear()
                continue # skip variants without a MeasureSet ID
            measureset_id = measureset[0].attrib['ID']
            mutant_allele = 'ALT' # default is that each entry refers to the alternate allele
            attributes = elem.findall('.//Attribute')
            for attribute in attributes:
                attribute_type = attribute.attrib.get('Type')
                if attribute_type is not None and "HGVS" in attribute_type and "protein" not in attribute_type: # if this is an HGVS cDNA, _not_ protein, annotation:
                    if attribute.text is not None and "=" in attribute.text: # and if there is an equals sign in the text, then
                        mutant_allele = 'REF' # that is their funny way of saying this assertion refers to the reference allele
            # find all the Citation nodes, and get the PMIDs out of them
            citations = elem.findall('.//Citation')
            pmids = []
            for citation in citations:
                pmids += [id_node.text for id_node in citation.findall('.//ID') if id_node.attrib.get('Source')=='PubMed']
            # now find the Comment nodes, regex your way through the comments and extract anything that appears to be a PMID
            comments = elem.findall('.//Comment')
            comment_pmids = []
            for comment in comments:
                mentions_pubmed = re.search(mentions_pubmed_regex,comment.text)
                if mentions_pubmed is not None and mentions_pubmed.group(1) is not None:
                    remaining_text = mentions_pubmed.group(1)
                    while True:
                        pubmed_id_extraction = re.search(extract_pubmed_id_regex,remaining_text)
                        if pubmed_id_extraction is None:
                            break
                        elif pubmed_id_extraction.group(1) is not None:
                            comment_pmids += [pubmed_id_extraction.group(1)]
                            if pubmed_id_extraction.group(2) is not None:
                                remaining_text = pubmed_id_extraction.group(2)
            all_pmids = list(set(pmids + comment_pmids))
            # now find any/all submitters
            submitters = []
            submitter_nodes = elem.findall('.//ClinVarSubmissionID')
            for submitter_node in submitter_nodes:
                if submitter_node.attrib is not None and submitter_node.attrib.has_key('submitter'):
                    submitters.append(submitter_node.attrib['submitter'])
            # now find the disease(s) this variant is associated with
            traitsets = elem.findall('.//TraitSet')
            all_traits = []
            for traitset in traitsets:
                trait_type = ''
                trait_values = []
                if traitset.attrib is not None:
                    trait_type = str(traitset.attrib.get('Type'))
                    disease_name_nodes = traitset.findall('.//Name/ElementValue')
                    for disease_name_node in disease_name_nodes:
                        if disease_name_node.attrib is not None:
                            if disease_name_node.attrib.get('Type') == 'Preferred':
                                trait_values.append(disease_name_node.text)
                all_traits += trait_values
            # now we're done traversing that one clinvar set. print out a cartesian product of accessions and pmids
            dest.write(('\t'.join( [chrom, pos, ref, alt, mutant_allele, measureset_id, ';'.join(submitters), ';'.join(all_traits), ','.join(all_pmids)] ) + '\n').encode('utf-8'))
            counter += 1
            if counter % 100 == 0:
                dest.flush()
            if verbose:
                sys.stderr.write("{0} entries completed, {1}, {2} total \r".format(
                    counter, ', '.join('%s skipped due to %s' % (v, k) for k, v in skipped_counter.items()),
                    counter + sum(skipped_counter.values())))
                sys.stderr.flush()
            elem.clear()

def get_handle(path):
    if path[-3:] == '.gz':
        handle = gzip.open(path)
    else:
        handle = open(path)
    return (handle)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract PMIDs from the ClinVar XML dump')
    parser.add_argument('-x','--xml', dest='xml_path',
                       type=str, help='Path to the ClinVar XML dump')
    parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),
                       default=sys.stdout)
    args = parser.parse_args()
    parse_clinvar_tree(get_handle(args.xml_path),dest=args.out)
