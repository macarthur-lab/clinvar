### ClinVar

#### In 1 sentence

This branch modifies the [clinvar](https://github.com/macarthur-lab/clinvar) master branch to convert latest ClinVar xml file into a tab-delimited file.

#### New Features

1. Extract the the variant info and interpretation info e.g. clinical significance and review status from ClinVar xml file for the measure set 
2. Record the measure set type in order to differentiate the haplotype with more than one alleles.


#### Usage
```
cd ./src
pip install --user --upgrade -r requirements.txt
python master_parse_clinvar_xml.py -R hg19.fa
```

#### Results

The main output files are:
* [clinvar.tsv.gz](clinvar.tsv.gz)  


