'''
Created on 16 Dec 2016

@author: xzhang13
'''
'''
Created on 23 Nov 2016

@author: xzhang13
'''
#modify the clinvar package code in order to parse xml to table with with more fields.


import re
import sys
import gzip
import argparse
from collections import defaultdict
import xml.etree.ElementTree as ET

# then sort it: cat clinvar_table.tsv | head -1 > clinvar_table_sorted.tsv; cat clinvar_table.tsv | tail -n +2 | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv


#Reference: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/clinvar_submission.xsd
#Reference on the XML info: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README

mentions_pubmed_regex = '(?:PubMed|PMID)(.*)' # group(1) will be all the text after the word PubMed or PMID
extract_pubmed_id_regex = '[^0-9]+([0-9]+)[^0-9](.*)' # group(1) will be the first PubMed ID, group(2) will be all remaining text

def replace_semicolons(s, replace_with=":"):
    return s.replace(";", replace_with)

def remove_newlines_and_tabs(s):
    return re.sub("[\t\n\r]", " ", s)

def parse_clinvar_tree(handle,dest=sys.stdout,error=sys.stderr,verbose=True,mode='collapsed'):
    # print a header row
    #header = [
    #    'chrom', 'pos', 'ref', 'alt', 'mut', 'alleleid', 'name','all_submitters', 'all_traits', 'all_pmids',
    #    'inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'xrefs'
    #]
    header = [
        'chrom', 'pos', 'ref', 'alt', 'mut', 'Measureset_type','Measureset_id','number_alleles','Measure_type','alleleid', 'name','gene_symbol',
        'molecular_consequence','variant_type','variant_type_change','all_conditions_name'
        ,'clinical_significance','review_status','all_submitters'
        , 'all_pmids','inheritance_modes', 'age_of_onset', 'prevalence', 
        'disease_mechanism', 'xrefs'
    ]
    dest.write(('\t'.join(header) + '\n').encode('utf-8'))
    counter = 0
    skipped_counter = defaultdict(int)
    for event, elem in ET.iterparse(handle):
        if elem.tag != 'ClinVarSet' or event != 'end':
            continue
        
        #initialize all the fields

        #Xiaolei
        #since the same measureset can have two measure alleles, which can not be interpretated independently thus reported together as haplotype or compound heterzygous. 
        #for this case, record the variation type in Measureset_type, and if the number_alleles is larger than one, aggregate all the other fields.  
        
        current_row = {}
        current_row['Measureset_type']=''
        current_row['Measureset_id']=0
        current_row['number_alleles']=0
        current_row['Measure_type']={}
        current_row['alleleid']={}
        
        measureset = elem.findall(".//ReferenceClinVarAssertion/MeasureSet")
        if(len(measureset)>1):
            print "A submission has more than one measure set."+elem.find('./Title').text #which is not possible?
            elem.clear()
            continue
        elif(len(measureset)==0):
            print "A submission has no measure set type"+measureset.attrib.get('ID')
            elem.clear()
            continue
        else:
            measureset=measureset[0]
            current_row['Measureset_id']=measureset.attrib.get('ID')
            current_row['Measureset_type']=measureset.get('Type')
            current_row['number_alleles']=str(len(measureset.findall('.//Measure')))
        
        #only the ones with just one measure set can be recorded
        
        
        measure_types=[]
        alleleids=[]
        #Aggregate all the possible measure types
        for measure_type in measureset.findall('.//Measure'):
            measure_types.append(measure_type.attrib.get('Type'))
            alleleids.append(measure_type.attrib['ID'])
        
        current_row['Measure_type']=";".join(measure_types)
        current_row['alleleid']=";".join(alleleids)
        
      
        #Xiaolei: Even the submission is more than one alleles, still want to keep it for reference to seperate whether there are independent interpretation in other submissions. 
        #Xiaolei: If there are more than one alleles in one submission, in the column fields (Chrom Pos REF ALT MUT gene_symbol), one of allele (the first one) information would be recorded. 
        #Xiaolei: For the other fields, the entries are aggregated
            
        # find the GRCh37 VCF representation
        grch37_location = None
        for sequence_location in elem.findall(".//SequenceLocation"):
            if sequence_location.attrib.get('Assembly') == 'GRCh37':
                if all(sequence_location.attrib.get(key) is not None for key in ('Chr', 'start', 'referenceAllele','alternateAllele')):
                    grch37_location = sequence_location
                    #print "need to break since the sequence"+current_row['Measureset_id'] 
                    break;
        #break after finding the first non-empty GRCh37 location
                
        if grch37_location is None:
            skipped_counter['missing SequenceLocation'] += 1
            elem.clear()
            continue # don't bother with variants that don't have a VCF location
        
        #find the allele ID (//Measure/@ID)
        measure = elem.findall('.//Measure')
        if measure is None:
            skipped_counter['missing Measure'] += 1
            elem.clear()
            continue # skip variants without a Measure ID

        
        current_row['chrom'] = grch37_location.attrib['Chr']
        current_row['pos'] = grch37_location.attrib['start']
        current_row['ref'] = grch37_location.attrib['referenceAllele']
        current_row['alt'] = grch37_location.attrib['alternateAllele']


        # iterate over attributes in the MeasureSet
        current_row['mut'] = 'ALT' # default is that each entry refers to the alternate allele
        for attribute in elem.findall('.//Attribute'):
            attribute_type = attribute.attrib.get('Type')
            if attribute_type is not None and "HGVS" in attribute_type and "protein" not in attribute_type: # if this is an HGVS cDNA, _not_ protein, annotation:
                if attribute.text is not None and "=" in attribute.text: # and if there is an equals sign in the text, then
                    current_row['mut'] = 'REF' # that is their funny way of saying this assertion refers to the reference allele

        # find the name HGVS representation http://varnomen.hgvs.org
        name = elem.findall('./ReferenceClinVarAssertion/MeasureSet/Name')
        if name is None:
            skipped_counter['missing variant name'] += 1
            elem.clear()
            continue # skip variants without a name
        for name_id in name:
            for ElementValue in name_id.findall('.//ElementValue'):
                if(ElementValue.attrib.get('Type')=="Preferred"):
                    current_row['name']=ElementValue.text;
                    break
            
        #find the gene symbol 
        current_row['gene_symbol']=''
        genesymbol = elem.findall('.//Symbol')
        if not genesymbol:
            skipped_counter['missing gene symbol'] += 1
            elem.clear()
            continue # skip variants without a gene symbol
        if(genesymbol[0].find('ElementValue').attrib.get('Type')=='Preferred'):
            current_row['gene_symbol']=genesymbol[0].find('.//ElementValue').text;

            
            
        attributeset=elem.findall('./ReferenceClinVarAssertion/MeasureSet/Measure/AttributeSet')
        current_row['molecular_consequence']=set()
        exist_variant_type=False;
        current_row['variant_type']=''
        current_row['variant_type_change']=''
        
        for attribute_node in attributeset: 
            attribute_type=attribute_node.find('./Attribute').attrib.get('Type')
            attribute_value=attribute_node.find('./Attribute').text;
            
            #Optional field
            #aggregate all molecular consequence
            if (attribute_type=='MolecularConsequence'):
                for xref in attribute_node.findall('.//XRef'):
                    if(xref.attrib.get('DB')=="RefSeq"):
                        #print xref.attrib.get('ID'), attribute_value
                        current_row['molecular_consequence'].add(":".join([xref.attrib.get('ID'),attribute_value]))
            
            
            #TODO:Check if it is mutually exclusive: protein vs nucleotide. Currently assumes that it is. 
            #currently just record the ones annotated as "nucleotide change" or "protein change"
            #Optional field 
            #find the variant_type: protein or nucleotide and variant_type_change: the change content
            if not exist_variant_type:
                exist_variant_type_nucleotide=(attribute_type=='nucleotide change')
                exist_variant_type_protein=(attribute_type=='ProteinChange1LetterCode' or attribute_type=='ProteinChange3LetterCode')
                if exist_variant_type_nucleotide is True or exist_variant_type_protein is True:
                    exist_variant_type=True
                    if exist_variant_type_nucleotide is True:
                        current_row['variant_type']='nucleotide change'
                        current_row['variant_type_change']=attribute_value
                    elif exist_variant_type_protein is True:
                        current_row['variant_type']='protein change'
                        current_row['variant_type_change'] = attribute_value
            
        
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
        #if it doesn't exist,skip the record
        if clinical_significance is None:
            skipped_counter['missing clicinical significance description'] += 1
            #elem.clear()
            continue # skip variants without a gene symbol
        if clinical_significance.find('.//ReviewStatus') is not None:
            current_row['review_status']=clinical_significance.find('.//ReviewStatus').text;
        if clinical_significance.find('.//Description') is not None:
            current_row['clinical_significance']=clinical_significance.find('.//Description').text
        
        # init new fields
        for list_column in ('inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'xrefs'):
            current_row[list_column] = set()

        # now find the disease(s) this variant is associated with
        current_row['all_conditions_name'] = []
        for traitset in elem.findall('.//TraitSet'):
            disease_name_nodes = traitset.findall('.//Name/ElementValue')
            trait_values = []
            for disease_name_node in disease_name_nodes:
                if disease_name_node.attrib is not None and disease_name_node.attrib.get('Type') == 'Preferred':
                    trait_values.append(disease_name_node.text)
            current_row['all_conditions_name'] += trait_values
            
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

        # done parsing the xml for this one clinvar set.
        elem.clear()

        # convert collection to string for the following fields
        for column_name in ('molecular_consequence','all_conditions_name', 'inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'xrefs'):
            column_value = current_row[column_name] if type(current_row[column_name]) == list else sorted(current_row[column_name])  # sort columns of type 'set' to get deterministic order
            current_row[column_name] = remove_newlines_and_tabs(';'.join(map(replace_semicolons, column_value)))

        # write out the current_row
#       :
        dest.write(('\t'.join([current_row[column] for column in header]) + '\n').encode('utf-8'))
#        except:
#            print current_row['alleleid']+'[alleleid] goes wrong'
            

        counter += 1
        if counter % 100 == 0:
            dest.flush()
        if verbose:
            error.write("{0} entries completed, {1}, {2} total \r".format(
                counter, ', '.join('%s skipped due to %s' % (v, k) for k, v in skipped_counter.items()),
                counter + sum(skipped_counter.values())))
            error.flush()
    error.write("Done\n")


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
    parser.add_argument('-e','--error', nargs='?', type=argparse.FileType('w'), default=sys.stderr)
    args = parser.parse_args()
    parse_clinvar_tree(get_handle(args.xml_path),dest=args.out,error=args.error)
