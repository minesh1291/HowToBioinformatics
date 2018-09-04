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

		
##
		
		
		
		
		
