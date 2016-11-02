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

def replace_semicolons(s, replace_with=":"):
    return s.replace(";", replace_with)

def remove_newlines_and_tabs(s):
    return re.sub("[\t\n\r]", " ", s)

def parse_clinvar_tree(handle,dest=sys.stdout,verbose=True,mode='collapsed'):
    # print a header row
    header = [
        'chrom', 'pos', 'ref', 'alt', 'mut', 'measureset_id', 'all_submitters', 'all_traits', 'all_pmids',
        'inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'xrefs'
    ]
    dest.write(('\t'.join(header) + '\n').encode('utf-8'))
    counter = 0
    skipped_counter = defaultdict(int)
    for event, elem in ET.iterparse(handle):
        if elem.tag != 'ClinVarSet' or event != 'end':
            continue

        # find the GRCh37 VCF representation
        grch37_location = None
        for sequence_location in elem.findall(".//SequenceLocation"):
            if sequence_location.attrib.get('Assembly') == 'GRCh37':
                if all(sequence_location.attrib.get(key) is not None for key in ('Chr', 'start', 'referenceAllele','alternateAllele')):
                    grch37_location = sequence_location

        if grch37_location is None:
            skipped_counter['missing SequenceLocation'] += 1
            elem.clear()
            continue # don't bother with variants that don't have a VCF location

        measuresets = elem.findall('.//MeasureSet')
        if measuresets is None:
            skipped_counter['missing MeasureSet'] += 1
            elem.clear()
            continue # skip variants without a MeasureSet ID

        current_row = {}
        current_row['chrom'] = grch37_location.attrib['Chr']
        current_row['pos'] = grch37_location.attrib['start']
        current_row['ref'] = grch37_location.attrib['referenceAllele']
        current_row['alt'] = grch37_location.attrib['alternateAllele']
        current_row['measureset_id'] = measuresets[0].attrib['ID']

        # iterate over attributes in the MeasureSet
        current_row['mut'] = 'ALT' # default is that each entry refers to the alternate allele
        for attribute in elem.findall('.//Attribute'):
            attribute_type = attribute.attrib.get('Type')
            if attribute_type is not None and "HGVS" in attribute_type and "protein" not in attribute_type: # if this is an HGVS cDNA, _not_ protein, annotation:
                if attribute.text is not None and "=" in attribute.text: # and if there is an equals sign in the text, then
                    current_row['mut'] = 'REF' # that is their funny way of saying this assertion refers to the reference allele

        # init list fields

        # find all the Citation nodes, and get the PMIDs out of them
        pmids = []
        for citation in elem.findall('.//Citation'):
            pmids += [id_node.text for id_node in citation.findall('.//ID') if id_node.attrib.get('Source') == 'PubMed']

        # now find the Comment nodes, regex your way through the comments and extract anything that appears to be a PMID
        comment_pmids = []
        for comment in elem.findall('.//Comment'):
            mentions_pubmed = re.search(mentions_pubmed_regex,comment.text)
            if mentions_pubmed is not None and mentions_pubmed.group(1) is not None:
                remaining_text = mentions_pubmed.group(1)
                while True:
                    pubmed_id_extraction = re.search(extract_pubmed_id_regex,remaining_text)
                    if pubmed_id_extraction is None:
                        break
                    elif pubmed_id_extraction.group(1) is not None:
                        comment_pmids.append( pubmed_id_extraction.group(1) )
                        if pubmed_id_extraction.group(2) is not None:
                            remaining_text = pubmed_id_extraction.group(2)

        current_row['all_pmids'] = ','.join(sorted(set(pmids + comment_pmids)))

        # now find any/all submitters
        current_row['all_submitters'] = ';'.join([
            submitter_node.attrib['submitter'].replace(';', ',')
            for submitter_node in elem.findall('.//ClinVarSubmissionID')
            if submitter_node.attrib is not None and submitter_node.attrib.has_key('submitter')
        ])

        # init new fields
        for list_column in ('inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'xrefs'):
            current_row[list_column] = set()

        # now find the disease(s) this variant is associated with
        current_row['all_traits'] = []
        for traitset in elem.findall('.//TraitSet'):
            disease_name_nodes = traitset.findall('.//Name/ElementValue')
            trait_values = []
            for disease_name_node in disease_name_nodes:
                if disease_name_node.attrib is not None and disease_name_node.attrib.get('Type') == 'Preferred':
                    trait_values.append(disease_name_node.text)
            current_row['all_traits'] += trait_values

            for attribute_node in traitset.findall('.//AttributeSet/Attribute'):
                attribute_type = attribute_node.attrib.get('Type')
                if attribute_type in {'ModeOfInheritance', 'age of onset', 'prevalence', 'disease mechanism'}:
                    column_name = 'inheritance_modes' if attribute_type == 'ModeOfInheritance' else attribute_type.replace(' ', '_')
                    column_value = attribute_node.text.strip()
                    if column_value:
                        current_row[column_name].add(column_value)

            for xref_node in traitset.findall('.//XRef'):
                xref_db = xref_node.attrib.get('DB')
                xref_id = xref_node.attrib.get('ID')
                current_row['xrefs'].add("%s:%s" % (xref_db, xref_id))

        # done parsing the xml for this one clinvar set.
        elem.clear()

        # convert collection to string for the following fields
        for column_name in ('all_traits', 'inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'xrefs'):
            column_value = current_row[column_name] if type(current_row[column_name]) == list else sorted(current_row[column_name])  # sort columns of type 'set' to get deterministic order
            current_row[column_name] = remove_newlines_and_tabs(';'.join(map(replace_semicolons, column_value)))

        # write out the current_row
        dest.write(('\t'.join([current_row[column] for column in header]) + '\n').encode('utf-8'))
        counter += 1
        if counter % 100 == 0:
            dest.flush()
        if verbose:
            sys.stderr.write("{0} entries completed, {1}, {2} total \r".format(
                counter, ', '.join('%s skipped due to %s' % (v, k) for k, v in skipped_counter.items()),
                counter + sum(skipped_counter.values())))
            sys.stderr.flush()
    sys.stderr.write("Done\n")


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
