rm(list=ls())
library(tidyverse)
expMatrix <- read.table("FastQ/GSE102349_featurecounts.txt",comment.char = "#",
                        row.names = 1, header = TRUE, sep="\t",check.names = F)
expMatrix[1:5,1:6]
# See the read count for the first five genes
exprSet <- expMatrix[,-c(1:5)]
exprSet[1:5,1:6]
colnames(exprSet) <- gsub("(.*?)\\/(.*?)\\_sorted.bam","\\2\\",colnames(exprSet) ) 
exprSet_out <- exprSet %>% as.data.frame()%>%  rownames_to_column("gene_id")
exprSet_out[1:5,1:6]
write.table(exprSet_out,file = "01rawData/GSE102349_featurecounts.txt",row.names = F,quote = F,sep = "\t")
# The count value is converted to TPM (Transcripts Per Million)  
x <- exprSet / expMatrix$Length
exprSet_tpm <- t( t(x) / colSums(x) ) * 1e6 
exprSet_tpm[1:3,1:5]
exprSet_tpm_out <- exprSet_tpm %>% as.data.frame()%>%  rownames_to_column("gene_id")
exprSet_tpm_out[1:3,1:5]
write.table(exprSet_tpm_out,file = "01rawData/GSE102349_featurecounts_TPM.txt",row.names = F,quote = F,sep = "\t")

# ensemble id is converted to symbol id 
library("rtracklayer")
gtf1 = import('D:/TCGA/Homo_sapiens.GRCh38.100.chr.gtf') 
gtf_df <- as.data.frame(gtf1)
gtf_df[1:9,1:3]
gtf_df <- gtf_df[gtf_df$gene_biotype == "protein_coding",c("gene_id","gene_biotype","gene_name")]
# gtf_df$"gene_id" <- substr(gtf_df$"gene_id", 1,15)
gtf_df <- unique(gtf_df)

exprSet_symbol <- exprSet_out %>%
  dplyr::inner_join(gtf_df,by="gene_id") %>%
  # tidyr::unite(gene_id,gene_id,gene_type,gene_name,sep = "|")%>%
  # Get rid of redundant information
  dplyr::select(-c(gene_id,gene_biotype))%>%
  # realignment
  dplyr::select(gene_name,everything()) %>%
  # The expression levels of duplicate genes were averaged
  dplyr::group_by(gene_name) %>%
  dplyr::summarise_all(mean)
exprSet_symbol[1:9,1:3]
write.table(exprSet_symbol,file = "01rawData/GSE102349_counts_symbol.txt",row.names = F,quote = F,sep = "\t")
save(exprSet_symbol,file = "01rawData/GSE102349_counts_symbol.Rdata")

# TPM symbol
exprSet_TPM_symbol <- exprSet_tpm_out %>%
  dplyr::inner_join(gtf_df,by="gene_id") %>%
  # tidyr::unite(gene_id,gene_id,gene_type,gene_name,sep = "|")%>%
  # Get rid of redundant information
  dplyr::select(-c(gene_id,gene_biotype))%>%
  #realignment
  dplyr::select(gene_name,everything()) %>%
  # The expression levels of duplicate genes were averaged
  dplyr::group_by(gene_name) %>%
  dplyr::summarise_all(mean)
exprSet_TPM_symbol[1:9,1:3]
write.table(exprSet_TPM_symbol,file = "01rawData/GSE102349_TPM_symbol.txt",row.names = F,quote = F,sep = "\t")
save(exprSet_TPM_symbol,file = "01rawData/GSE102349_TPM_symbol.Rdata")

