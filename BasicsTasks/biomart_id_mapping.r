# for getting ID of conventional transcript for given gene symbol

library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

filters = listFilters(ensembl)
filters[1:5,]

attributes = listAttributes(ensembl)
attributes[1:5,]

grep(pattern = "hsa.*cano",ignore.case = T,x=attributes$name,value = T) # hsapiens_paralog_canonical_transcript_protein
grep(pattern = "gene_id",ignore.case = T,x=attributes$name,value = T) #ensembl_gene_id
grep(pattern = ".*ensembl.*",ignore.case = T,x=attributes$name,value = T)

gs = c("PTPRC","IL2")

ipro = getBM(attributes=c("refseq_mrna","interpro","interpro_description"), 
						 filters="hgnc_symbol",
						 values=gs, 
						 mart=ensembl)

cononical_TrPr = getBM(attributes=c("hsapiens_paralog_canonical_transcript_protein"), 
											 filters="hgnc_symbol",
											 values=gs, 
											 mart=ensembl)
	#ENSP00000411355

gene2pep = getBM(attributes=c("hgnc_symbol","ensembl_gene_id","ensembl_peptide_id"), 
											 filters="hgnc_symbol",
											 values=gs, 
											 mart=ensembl)

gene2pep[gene2pep$ensembl_peptide_id=="ENSP00000411355",] 


# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownCanonical.txt.gz

# https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/  # good practice
