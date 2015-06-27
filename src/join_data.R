#!/usr/bin/env Rscript

options(stringsAsFactors=F)

# load what we've extracted from the XML so far
xml_extract = read.table('clinvar_table_dedup_normalized.tsv',sep='\t',header=T)

# load the tab-delimited summary
txt_download = read.table('variant_summary.txt.gz',sep='\t',comment.char='',quote='',header=T)

# subset the tab-delimited summary to desired rows and cols
colnames(txt_download) = gsub('\\.','_',tolower(colnames(txt_download)))
desired_columns = c('variantid','genesymbol','clinicalsignificance','reviewstatus','hgvs_c__','hgvs_p__')
txt_extract = subset(txt_download, assembly == 'GRCh37', select=desired_columns)
colnames(txt_extract) = c('measureset_id','symbol','clinical_significance','review_status','hgvs_c','hgvs_p')

# join on measureset_id
combined = merge(xml_extract, txt_extract, by='measureset_id')

# re-order the columns
combined = combined[,c('chrom','pos','ref','alt','mut','measureset_id','symbol','clinical_significance','review_status','hgvs_c','hgvs_p','all_submitters','all_pmids')]

# add some layers of interpretation on top of this
# note: we are trying to get the "overall" interpretation that is displayed in the upper right of the clinvar web pages but
# it is not in any of the available FTP downloads, so this is a stopgap
combined$pathogenic = as.integer(grepl('athogenic',combined$clinical_significance)) # 1 if at least one submission says path or likely path, 0 otherwise
combined$conflicted = as.integer(grepl('athogenic',combined$clinical_significance) & grepl('enign',combined$clinical_significance)) # 1 if at least one submission each of [likely] benign and [likely] pathogenic

write.table(combined,'clinvar_combined.tsv',sep='\t',row.names=F,col.names=T,quote=F)


