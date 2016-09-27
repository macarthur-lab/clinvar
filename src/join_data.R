#!/usr/bin/env Rscript

options(stringsAsFactors=F)
options(warn=2) 
options(error = quote({
  dump.frames(to.file=T, dumpto='last.dump')
  load('last.dump.rda')
  print(last.dump)
  q()
}))

# load what we've extracted from the XML so far
xml_extract = read.table('clinvar_table_normalized.tsv',sep='\t',comment.char='',quote='',header=T)

# load the tab-delimited summary
txt_download = read.table('variant_summary.txt.gz',sep='\t',comment.char='',quote='',header=T,fileEncoding="UTF-16LE")

# subset the tab-delimited summary to desired rows and cols
colnames(txt_download) = gsub('\\.','_',tolower(colnames(txt_download)))

desired_columns = c('variantid','genesymbol','clinicalsignificance','reviewstatus','hgvs_c__','hgvs_p__', 'origin')
txt_extract = subset(txt_download, assembly == 'GRCh37', select=desired_columns)
colnames(txt_extract) = c('measureset_id','symbol','clinical_significance','review_status','hgvs_c','hgvs_p','origin')


# join on measureset_id
combined = merge(xml_extract, txt_extract, by='measureset_id')

# re-order the columns
combined = combined[,c('chrom','pos','ref','alt','mut','measureset_id','symbol','clinical_significance','review_status','hgvs_c','hgvs_p','all_submitters','all_traits','all_pmids', 'inheritance_modes', 'age_of_onset', 'prevalence', 'disease_mechanism', 'origin', 'xrefs')]

# lookup table based on http://www.ncbi.nlm.nih.gov/clinvar/docs/details/
gold_stars_table = list(
	'no assertion provided' = 0,
	'no assertion for the individual variant' = 0,
	'no assertion criteria provided' = 0,
	'criteria provided, single submitter' = 1,
	'criteria provided, conflicting interpretations' = 1, 
	'criteria provided, multiple submitters, no conflicts' = 2, 
	'reviewed by expert panel' = 3,
	'practice guideline' = 4
)

# add some layers of interpretation on top of this
# note: we are trying to get the "overall" interpretation that is displayed in the upper right of the clinvar web pages but
# it is not in any of the available FTP downloads, so this is a stopgap
combined$gold_stars = sapply(combined$review_status, function(k) { gold_stars_table[[k]] })
# pathogenic = 1 if at least one submission says path or likely path, 0 otherwise
combined$pathogenic = as.integer(grepl('athogenic',combined$clinical_significance))
# conflicted = 1 if at least one submission each of [likely] benign and [likely] pathogenic
combined$conflicted = as.integer(grepl('athogenic',combined$clinical_significance) & grepl('enign',combined$clinical_significance))
# benign = 1 if at least one submission says benign or likely benign, 0 otherwise
combined$benign = as.integer(grepl('enign',combined$clinical_significance))

write.table(combined,'clinvar_combined.tsv',sep='\t',row.names=F,col.names=T,quote=F)
