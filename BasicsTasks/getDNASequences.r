# getting DNA sequences by GeneSymbols

library("BiocInstaller")
biocLite("biomaRt")
library(biomaRt)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

geneSymbols = c("PTPRC")

DNAseqs = getSequence(id = geneSymbols, 
						type='hgnc_symbol',
						seqType="gene_exon_intron",
						mart=ensembl)
