#!/usr/bin/env python

import re
import sys
import gzip
import argparse
from collections import defaultdict
import xml.etree.ElementTree as ET


# then sort it: cat clinvar_table.tsv | head -1 > clinvar_table_sorted.tsv; cat clinvar_table.tsv | tail -n +2 | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv
# Reference on clinvar XML tag: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/clinvar_submission.xsd
# Reference on clinvar XML tag: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README

mentions_pubmed_regex = '(?:PubMed|PMID)(.*)' # group(1) will be all the text after the word PubMed or PMID
extract_pubmed_id_regex = '[^0-9]+([0-9]+)[^0-9](.*)' # group(1) will be the first PubMed ID, group(2) will be all remaining text

def replace_semicolons(s, replace_with=":"):
    return s.replace(";", replace_with)

def remove_newlines_and_tabs(s):
    return re.sub("[\t\n\r]", " ", s)

def parse_clinvar_tree(handle, dest=sys.stdout, multi=None, verbose=True, genome_build='GRCh37'):
    """Parse clinvar XML
    Args:
        handle: Open input file handle for reading the XML data
        dest: Open output file handle or stream for simple variants
        multi: Open output file handle or stream for complex non-single-variant clinvar records
            (eg. compound het, haplotypes, etc.)
        verbose: Whether to write extra stats to stderr
        genome_build: Either 'GRCh37' or 'GRCh38'
    """

    #measureset -> rcv (one to many) 
    header = [
        'chrom', 'pos', 'ref', 'alt', 'measureset_type','measureset_id','rcv',
        'allele_id','symbol',
        'hgvs_c','hgvs_p','molecular_consequence',
        'clinical_significance','clinical_significance_ordered','review_status','review_status_ordered','all_submitters','all_traits',
        'all_pmids','inheritance_modes', 'age_of_onset', 'prevalence', 
        'disease_mechanism', 'origin','xrefs'
    ]
    dest.write(('\t'.join(header) + '\n').encode('utf-8'))
    if multi is not None:
        multi.write(('\t'.join(header) + '\n').encode('utf-8'))

    scounter = 0
    mcounter = 0
    skipped_counter = defaultdict(int)
    for event, elem in ET.iterparse(handle):
        if elem.tag != 'ClinVarSet' or event != 'end':
            continue
        
        #initialize all the fields
        current_row = {}
        current_row['rcv']=''
        current_row['measureset_type']=''
        current_row['measureset_id']=''
        current_row['allele_id']=''
        
        rcv = elem.find('./ReferenceClinVarAssertion/ClinVarAccession')
        if rcv.attrib.get('Type')!='RCV':
            print("Error, not RCV record")
            break
        else:
            current_row['rcv']=rcv.attrib.get('Acc')
        
        ReferenceClinVarAssertion=elem.findall(".//ReferenceClinVarAssertion")
        measureset = ReferenceClinVarAssertion[0].findall(".//MeasureSet")
        
        #only the ones with just one measure set can be recorded
        if len(measureset) > 1:
            print("A submission has more than one measure set."+elem.find('./Title').text)
            elem.clear()
            continue
        elif len(measureset) == 0:
            print("A submission has no measure set type"+measureset.attrib.get('ID'))
            elem.clear()
            continue
        
        measureset=measureset[0]

        measure=measureset.findall('.//Measure')
        
        current_row['measureset_id']=measureset.attrib.get('ID')
        current_row['measureset_type']=measureset.get('Type')
        
        
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

        #find the clincial significance and review status reported in RCV(aggregated from SCV)
        current_row['clinical_significance']=[]
        current_row['review_status']=[]
        
        clinical_significance=elem.find('.//ReferenceClinVarAssertion/ClinicalSignificance') 
        if clinical_significance.find('.//ReviewStatus') is not None:
            current_row['review_status']=clinical_significance.find('.//ReviewStatus').text;
        if clinical_significance.find('.//Description') is not None:
            current_row['clinical_significance']=clinical_significance.find('.//Description').text

        #match the order of the submitter list - edit 2/22/17
        current_row['review_status_ordered'] = ';'.join([
            x.text for x in elem.findall('.//ClinVarAssertion/ClinicalSignificance/ReviewStatus') if x is not None
        ])
        current_row['clinical_significance_ordered'] = ';'.join([
            x.text for x in elem.findall('.//ClinVarAssertion/ClinicalSignificance/Description') if x is not None
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
                        
        #put all the cross references one column, it may contains NCBI gene ID, conditions ID in disease databases. 
            for xref_node in traitset.findall('.//XRef'):
                xref_db = xref_node.attrib.get('DB')
                xref_id = xref_node.attrib.get('ID')
                current_row['xrefs'].add("%s:%s" % (xref_db, xref_id))
        
        current_row['origin']=set()
        for origin in elem.findall('.//ReferenceClinVarAssertion/ObservedIn/Sample/Origin'):
            current_row['origin'].add(origin.text)
        
        for column_name in ('all_traits', 'inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'origin','xrefs'):
            column_value = current_row[column_name] if type(current_row[column_name]) == list else sorted(current_row[column_name])  # sort columns of type 'set' to get deterministic order
            current_row[column_name] = remove_newlines_and_tabs(';'.join(map(replace_semicolons, column_value)))

        for i in range(len(measure)):
            #find the allele ID (//Measure/@ID)
            current_row['allele_id']=measure[i].attrib.get('ID')
            # find the GRCh37 or GRCh38 VCF representation
            genomic_location = None
            for sequence_location in measure[i].findall(".//SequenceLocation"):
                if sequence_location.attrib.get('Assembly') == genome_build:
                    if all(sequence_location.attrib.get(key) is not None for key in ('Chr', 'start', 'referenceAllele','alternateAllele')):
                        genomic_location = sequence_location
                        break

            #break after finding the first non-empty GRCh37 or GRCh38 location
            if genomic_location is None:
                skipped_counter['missing SequenceLocation'] += 1
                elem.clear()
                continue # don't bother with variants that don't have a VCF location
                
            current_row['chrom'] = genomic_location.attrib['Chr']
            current_row['pos'] = genomic_location.attrib['start']
            current_row['ref'] = genomic_location.attrib['referenceAllele']
            current_row['alt'] = genomic_location.attrib['alternateAllele']
           
            #find the gene symbol 
            current_row['symbol']=''
            genesymbol = measure[i].findall('.//Symbol')
            if genesymbol is not None:
                for symbol in genesymbol:
                    if(symbol.find('ElementValue').attrib.get('Type')=='Preferred'):
                        current_row['symbol']=symbol.find('ElementValue').text;
                        break

            current_row['molecular_consequence']=set()
            current_row['hgvs_c']=''
            current_row['hgvs_p']=''

            attributeset=measure[i].findall('./AttributeSet')
            for attribute_node in attributeset: 
                attribute_type=attribute_node.find('./Attribute').attrib.get('Type')
                attribute_value=attribute_node.find('./Attribute').text;
            
                #find hgvs_c
                if(attribute_type=='HGVS, coding, RefSeq'):
                    current_row['hgvs_c']=attribute_value
            
                #find hgvs_p
                if(attribute_type=='HGVS, protein, RefSeq'):
                    current_row['hgvs_p']=attribute_value
            
                #aggregate all molecular consequences
                if (attribute_type=='MolecularConsequence'):
                    for xref in attribute_node.findall('.//XRef'):
                        if xref.attrib.get('DB')=="RefSeq":
                        #print xref.attrib.get('ID'), attribute_value
                            current_row['molecular_consequence'].add(":".join([xref.attrib.get('ID'),attribute_value]))
                
            column_name = 'molecular_consequence'
            column_value = current_row[column_name] if type(current_row[column_name]) == list else sorted(current_row[column_name])  # sort columns of type 'set' to get deterministic order
            current_row[column_name] = remove_newlines_and_tabs(';'.join(map(replace_semicolons, column_value)))
            
            if len(measure)==1:
                dest.write(('\t'.join([current_row[column] for column in header]) + '\n').encode('utf-8'))
                scounter += 1
            else:
                if multi is not None:
                    multi.write(('\t'.join([current_row[column] for column in header]) + '\n').encode('utf-8'))
                    mcounter += 1
            
            if scounter % 100 == 0:
                dest.flush()
            if mcounter % 100 == 0:
                if multi is not None:
                    multi.flush()
            
            counter = scounter + mcounter
            if verbose:
                sys.stderr.write("{0} entries completed, {1}, {2} total \r".format(
                    counter,
                    ', '.join('%s skipped due to %s' % (v, k) for k, v in skipped_counter.items()),
                    counter + sum(skipped_counter.values())
                ))
                sys.stderr.flush()
        
        # done parsing the xml for this one clinvar set.
        elem.clear()

    sys.stderr.write("Done\n")


def get_handle(path):
    if path[-3:] == '.gz':
        handle = gzip.open(path)
    else:
        handle = open(path)
    return handle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract PMIDs from the ClinVar XML dump')
    parser.add_argument('-g','--genome-build', choices=['GRCh37', 'GRCh38'],
                        help='Genome version (either GRCh37 or GRCh38)', required=True)
    parser.add_argument('-x','--xml', dest='xml_path',
                       type=str, help='Path to the ClinVar XML dump', required=True)
    parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-m', '--multi', help="Output file name for complex alleles")

    args = parser.parse_args()
    if args.multi is not None:
        f = open(args.multi, 'w')
        parse_clinvar_tree(get_handle(args.xml_path), dest=args.out, multi=f, genome_build=args.genome_build)
        f.close()
    else:
        parse_clinvar_tree(get_handle(args.xml_path), dest=args.out, genome_build=args.genome_build)
