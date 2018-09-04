# install.packages("bigrquery")
library(bigrquery)

## set billing and project
		Sys.setenv(BIGQUERY_TEST_PROJECT="tcga-data-analysis-")
		billing <- bq_test_project()
		set_service_token(service_token = "~/keys/tcga-data-analysis-.json")


## format query
		sql <- "SELECT year, month, day, weight_pounds FROM `publicdata.samples.natality`"
## fire query and download result
		tb <- bq_project_query(billing, sql,dryRun = F,quiet = F)
		tb_data <- bq_table_download(tb,max_results = 20)
		
		
## format query
		sql <- "select * FROM `isb-cgc:tcga_cohorts.__TABLES__` "
## fire query and download result
		tb <- bq_project_query(billing, sql,dryRun = F,quiet = F)
		tb_data <- bq_table_download(tb)
		tcga_cohorts <- as.data.frame(tb_data)
		save(tcga_cohorts,file = "tcga_cohorts.Rdata")		

# > tcga_cohorts
#    project_id   dataset_id table_id creation_time last_modified_time row_count size_bytes type
# 1     isb-cgc tcga_cohorts      ACC            NA                 NA       180       5760    1
# 2     isb-cgc tcga_cohorts     BLCA            NA                 NA       806      25792    1
# 3     isb-cgc tcga_cohorts     BRCA            NA                 NA      2236      71552    1
# 4     isb-cgc tcga_cohorts     CESC            NA                 NA       597      19104    1
		
## format query
		sql <- "SELECT
  gene_name,
  gene_type,
  Ensembl_gene_id_v
FROM
  `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression`
GROUP BY
  gene_name,
  gene_type,
  Ensembl_gene_id_v"
## fire query and download result
		tb <- bq_project_query(billing, query =  sql,dryRun = F,quiet = F)
		tb_data <- bq_table_download(tb)
		mRNA_genes <- as.data.frame(tb_data)
		# dim(mRNA_genes) # 60483
		# c("GLRX3","PTPRC") %in% mRNA_genes$gene_name
		# https://www.gencodegenes.org/gencode_biotypes.html for gene_type
		table(mRNA_genes$gene_type)
		save(mRNA_genes,file = "mRNA_genes.Rdata")		

#                              GeneType  Freq
# 20                     protein_coding 19814
# 18               processed_pseudogene 10304
# 10                            lincRNA  7656
# 2                           antisense  5565
# 12                              miRNA  4093
# 43             unprocessed_pseudogene  2574
# 13                           misc_RNA  2298
# 28                              snRNA  1896
# 30                                TEC  1045
		
##
		
		
		
		
		
