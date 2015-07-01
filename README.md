### ClinVar

#### In 1 sentence

This repo provides tools to convert ClinVar data into a tab-delimited flat file, and also provides that resulting tab-delimited flat file.

#### Motivation

[ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) is a public database hosted by NCBI for the purpose of collecting assertions as to genotype-phenotype pairings in the human genome. One common use case for ClinVar is as a catalogue of genetic variants that have been reported to cause Mendelian disease. In our work in the [MacArthur Lab](http://macarthurlab.org/), we have two major use cases for ClinVar:

1. To check whether candidate causal variants we find in Mendelian disease exomes have been previously reported as pathogenic.
2. To pair with [ExAC](http://exac.broadinstitute.org/) data to enable exome-wide analyses of reportedly pathogenic variants.

ClinVar makes its data available via [FTP](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/) in three formats: XML, TXT, and VCF. We found that none of these files were ideally suited for our purposes. The VCF only contains variants present in dbSNP; it is not a comprehensive catalogue of ClinVar variants. The TXT file lacks certain annotations such as PubMed IDs for related publications. The XML file is large and complex, with multiple entries for the same genomic variant, making it difficult to quickly look up a variant of interest. In addition, both the XML and TXT representations are not guaranteed to be unique on genomic coordinates, and also contain many genomic coordinates that have been parsed from HGVS notation, and therefore may be right-aligned (in contrast to left alignment, the standard for VCF) and may also be non-minimal (containing additional nucleotides of context to the left or right of a given variant).

#### Solution

To create a flat representation of ClinVar suited for our purposes, we took several steps, which are encapsulated in [src/master.bash](src/master.bash):

1. Download the latest XML and TXT dumps from ClinVar FTP.
2. Parse the XML file using [src/parse_clinvar_xml.py](src/parse_clinvar_xml.py) to extract fields of interest into a flat file.
3. Sort on genomic coordinates.
4. De-duplicate using [src/dedup_clinvar.py](src/dedup_clinvar.py), combining records that refer to the same genomic variant.
5. Normalize using [our Python implementation](https://github.com/ericminikel/minimal_representation/blob/master/normalize.py) of [vt normalize](http://genome.sph.umich.edu/wiki/Variant_Normalization) (see [[Tan 2015]]).
6. Join to some of the fields of interest from the TXT file using [src/join_data.R](src/join_data.R), and create some new fields&dagger;.
7. Sort and de-duplicate again (this removes dups arising from duplicate records in the TXT dump).

&dagger;Because a ClinVar record may contain multiple assertions of Clinical Significance, we defined two additional columns:

+ `pathogenic` is `1` if the variant has *ever* been asserted "Pathogenic" or "Likely pathogenic" by any submitter for any phenotype, and `0` otherwise
+ `conflicted` is `1` if the variant has *ever* been asserted "Pathogenic" or "Likely pathogenic" by any submitter for any phenotype, and has also been asserted "Benign" or "Likely benign" by any submitter for any phenotype, and `0` otherwise. Note that having one assertion of pathogenic and one of uncertain significance does *not* count as conflicted for this column. 

#### Results

The resulting output file is at [output/clinvar.tsv](output/clinvar.tsv).

#### Usage notes

Because ClinVar contains a great deal of data complexity, we made a deliberate decision to *not* attempt to capture all fields in our resulting file. We made an effort to capture a subset of fields that we believed would be most useful for genome-wide filtering, and also included `measureset_id` as a column to enable the user to look up additional details on the ClinVar website. For instance, the page for the variant with `measureset_id` 7105 is located at [ncbi.nlm.nih.gov/clinvar/variation/7105/](http://www.ncbi.nlm.nih.gov/clinvar/variation/7105/).

[Tan 2015]: http://www.ncbi.nlm.nih.gov/pubmed/25701572 "Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants. Bioinformatics. 2015 Jul 1;31(13):2202-4. doi: 10.1093/bioinformatics/btv112. Epub 2015 Feb 19. PubMed PMID: 25701572."