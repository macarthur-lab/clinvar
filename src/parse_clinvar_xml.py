#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

import re
import sys
import gzip
import argparse
import itertools
import xml.etree.ElementTree as ET

# to test: ./parse_clinvar_xml.py -x clinvar_test.xml
# to run for reals: bsub -q week -R rusage[mem=32] -oo cvxml.o -eo cvxml.e -J cvxml "./parse_clinvar_xml.py -x ClinVarFullRelease_2015-04.xml.gz -o clinvar_table.tsv"
# then sort it: cat clinvar_table.tsv | head -1 > clinvar_table_sorted.tsv; cat clinvar_table.tsv | tail -n +2 | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv

# the ClinVar XML dump is enormous - 113 MB gzipped, 2 GB gunzipped
# therefore, first subset to a reasonable number of entries to debug code:
# $ zcat ClinVarFullRelease_2015-04.xml.gz | grep -n -m 100 "/ClinVarSet"
# 40273 # is the line number of the end of the 100th entry, therefore:
# $ zcat ClinVarFullRelease_2015-04.xml.gz | head -40274 > clinvar_test.xml
# $ zcat ClinVarFullRelease_2015-04.xml.gz | tail -1 >> clinvar_test.xml

# also: use this to find lines where PMIDs are referenced:
# $ zcat ClinVarFullRelease_2015-04.xml.gz | egrep -m 2 -n "(PubMed|PMID)"

# pmid_regex = '(?:PubMed|PMID)[^0-9]*([0-9]+)[^0-9](.*)' # group(1) will be the first PMID, group(2) will be all remaining text
# pmid_regex_step_2 = '[^0-9]+?([0-9]+)[^0-9](.*)' # group(1) wi

mentions_pubmed_regex = '(?:PubMed|PMID)(.*)' # group(1) will be all the text after the word PubMed or PMID
extract_pubmed_id_regex = '[^0-9]+([0-9]+)[^0-9](.*)' # group(1) will be the first PubMed ID, group(2) will be all remaining text

revstat_ranking = {'not classified by submitter': 0, 
'classified by single submitter': 1, 
'classified by multiple submitters': 2, 
'reviewed by expert panel': 3,
'reviewed by professional society': 4,
}

# the full list of clinical significance Description values
clnsig_ranking = {'Uncertain significance': 0,
'not provided': 1,
'Benign': 2,
'Likely benign': 3,
'Likely pathogenic': 4,
'Pathogenic': 5,
'drug response': 6,
'histocompatibility': 7,
'other': 255}

# the subset of that list which we use for ranking among entries tied in ReviewStatus
# note though officially these are Sentence case http://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
# I find in practice they are Sentence case in some entries and all-lowercase in others
# so I have made them lowercase in this dictionary and then i cast everything to lowercase in the code below.
# also, Anne wants VUS to be 3.5 rather than 0 which it is above.
clnsig_ranking_path_spectrum = {
'benign': 2,
'likely benign': 3,
'uncertain significance': 3.5,
'likely pathogenic': 4,
'pathogenic': 5}

def parse_clinvar_tree(xml_tree,dest=sys.stdout,verbose=True,mode='collapsed'):
    # print a header row
    dest.write(('\t'.join( ['chrom', 'pos', 'ref', 'alt', 'accession', 'conflicted', 'mut', 'title', 'highest_revstat', 'worst_assertion', 'all_traits', 'all_submitters', 'all_clnsig', 'all_pmids'] ) + '\n').encode('utf-8'))
    root = xml_tree.getroot()
    clinvarsets = root.findall('ClinVarSet')
    counter = 0
    for clinvarset in clinvarsets:
        # find the GRCh37 VCF representation
        sequence_locations = clinvarset.findall('.//SequenceLocation')
        grch37 = None
        for sequence_location in sequence_locations:
            if sequence_location.attrib.get('Assembly') == 'GRCh37':
                if all(entry is not None for entry in [sequence_location.attrib.get(key) for key in ['referenceAllele','alternateAllele','start','Chr']]):
                    grch37 = sequence_location
        if grch37 is None:
            continue # don't bother with variants that don't have a VCF location
        else:
            chrom = grch37.attrib['Chr']
            pos = grch37.attrib['start']
            ref = grch37.attrib['referenceAllele']
            alt = grch37.attrib['alternateAllele']
        # get the human-readable title of the clinvar entry
        titles = clinvarset.findall('.//Title')
        if len(titles) == 0:
            title = ''
        else:
            title = titles[0].text
        # get all the accession numbers
        accession_numbers = [accession.attrib['Acc'] for accession in clinvarset.findall('.//ClinVarAccession')]
        accession_numbers = [x for x in accession_numbers if x[:3] == 'RCV'] # limit to RCV ones
        # get the clinical significance
        all_clnsig = []
        clnsigs = clinvarset.findall('.//ClinicalSignificance')
        max_revstat = 0 # highest review status (0 through 4)
        max_revstat_string = '' # string of highest review status (e.g. "classified by single submitter")
        max_description = 2 # highest pathogenicity assertion (2 through 5)
        max_description_string = '' # string, e.g. "Pathogenic"
        conflicted = False # is the entry conflicted? default is no.
        descriptions = [] # this will be a list of pathogenicity assertions, e.g. "likely pathogenic"
        for clnsig in clnsigs:
            # store the full text
            revstat_node = clnsig.find('ReviewStatus')
            if revstat_node is None:
                continue
            else:
                revstat  = revstat_node.text.lower() # convert to lowercase
            description_node = clnsig.find('Description')
            if description_node is None:
                continue
            else:
                description = description_node.text.lower() # convert to lowercase
            descriptions.append(description)
            all_clnsig.append(description + " (" + revstat + ")")
            # find the highest ranking entry
            if revstat_ranking.get(revstat, 0) >= max_revstat:
                max_revstat = revstat_ranking.get(revstat, 0)
                max_revstat_string = revstat
                if clnsig_ranking_path_spectrum.get(description, 0) >= max_description:
                    max_description = clnsig_ranking_path_spectrum.get(description, 0)
                    max_description_string = description
        if len(set(descriptions)) > 1:
            conflicted = True
        all_clnsig_text = ';'.join(all_clnsig)
        # figure out if the reference allele is the mutation
        # frustratingly, the _only_ way this is encoded in the ClinVar XML is that the
        # HGVS c. notation has an equals sign in it
        mutant_allele = 'ALT' # default is that each entry refers to the alternate allele
        attributes = clinvarset.findall('.//Attribute')
        for attribute in attributes:
            attribute_type = attribute.attrib.get('Type')
            if attribute_type is not None and "HGVS" in attribute_type and "protein" not in attribute_type: # if this is an HGVS cDNA, _not_ protein, annotation:
                if attribute.text is not None and "=" in attribute.text: # and if there is an equals sign in the text, then
                    mutant_allele = 'REF' # that is their funny way of saying this assertion refers to the reference allele
        # find all the Citation nodes, and get the PMIDs out of them
        citations = clinvarset.findall('.//Citation')
        pmids = []
        for citation in citations:
            pmids += [id_node.text for id_node in citation.findall('.//ID') if id_node.attrib.get('Source')=='PubMed']
        # now find the Comment nodes, regex your way through the comments and extract anything that appears to be a PMID
        comments = clinvarset.findall('.//Comment')
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
        # now find the disease(s) this variant is associated with
        traitsets = clinvarset.findall('.//TraitSet')
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
        all_traits_text = '; '.join(all_traits)
        # now find any/all submitters
        submitters = []
        submitter_nodes = clinvarset.findall('.//ClinVarSubmissionID')
        for submitter_node in submitter_nodes:
            if submitter_node.attrib is not None and submitter_node.attrib.has_key('submitter'):
                submitters.append(submitter_node.attrib['submitter'])
        all_submitters_text = '; '.join(submitters)
        # now we're done traversing that one clinvar set. print out a cartesian product of accessions and pmids
        if mode == 'relational':
            for element in itertools.product(accession_numbers,all_pmids):
                dest.write('\t'.join(element)+'\n')
        elif mode == 'collapsed':
            dest.write(('\t'.join( [chrom, pos, ref, alt, ','.join(accession_numbers), str(int(conflicted)), mutant_allele, title, max_revstat_string, max_description_string, all_traits_text, all_submitters_text, all_clnsig_text, ','.join(all_pmids)] ) + '\n').encode('utf-8'))
        counter += 1
        if counter % 100 == 0:
            dest.flush()
        if verbose:
            sys.stderr.write("\r{0} entries completed".format(counter))
            sys.stderr.flush()

def open_clinvar_tree(path_to_clinvar):
    if path_to_clinvar[-3:] == '.gz':
        handle = gzip.open(path_to_clinvar)
    else:
        handle = open(path_to_clinvar)
    tree = ET.parse(handle)
    return tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract PMIDs from the ClinVar XML dump')
    parser.add_argument('-x','--xml', dest='xml_path',
                       type=str, help='Path to the ClinVar XML dump')
    parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),
                       default=sys.stdout)
    args = parser.parse_args()
    tree = open_clinvar_tree(args.xml_path)
    parse_clinvar_tree(tree,dest=args.out)
