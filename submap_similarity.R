rm(list = ls())  
options(stringsAsFactors = F)
expFile_A="01rawData/GSE102349_TPM_symbol.txt"                                                 
pdFile_A = "01rawData/GSE102349_Cluster_res.txt"

# Read the expression result file and organize the data
testExpr_A <- read.table(expFile_A,sep="\t",header=T,row.names=1,check.names=F)
pd_A<- read.table(pdFile_A,sep="\t",header=T,check.names=F)
# Cluster1_A <- pd_A$id[pd_A$Cluster =="Cluster1"]
# Cluster2_A <- pd_A$id[pd_A$Cluster =="Cluster2"]

expFile_B="01rawData/GSE53819_NPC_exprSet.txt"                                                 
pdFile_B = "01rawData/GSE53819_Cluster_res.txt"

# Read the expression result file and organize the data
testExpr_B <- read.table(expFile_B,sep="\t",header=T,row.names=1,check.names=F)
pd_B <- read.table(pdFile_B,sep="\t",header=T,check.names=F)
Cluster1_B <- pd_B$id[pd_B$Cluster =="Cluster1"]
Cluster2_B <- pd_B$id[pd_B$Cluster =="Cluster2"]


 
# Custom functions are used to generate the data format required by submap
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# The data format needed to create the submap (expFile_B)
skcm.immunotherapy.logNC <- read.table(expFile_B,sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) 
skcm.immunotherapy.info <- read.table(pdFile_B,sep = "\t",header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$Cluster),]
skcm.immunotherapy.info$rank <- rep(c(1,2),times=as.character(table(skcm.immunotherapy.info$Cluster))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

# The data format needed to create the submap (expFile_A)
tmp <- read.table(expFile_A,sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1) 
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # Take the list of genes after the intersection of two data sets

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,skcm.immunotherapy.info$id]


# The file name that generates the output data
gct_file <- "GSE53819.for.SubMap.gct"
cls_file <- "GSE53819.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# Propose subtypes of samples and arrange them sequentially
samples.C1 <- pd_A$id[which(pd_A$Cluster == "Cluster1")]
samples.C2 <- pd_A$id[which(pd_A$Cluster == "Cluster2")]

sam_info <- data.frame("Cluster"=c(samples.C1,samples.C2),row.names = c(samples.C1,samples.C2))
sam_info$rank <- rep(c(1,2),times=c(length(samples.C1),length(samples.C2))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

# The file name that generates the output data
gct_file <- "GSE102349.for.SubMap.gct"
cls_file <- "GSE102349.for.SubMap.cls"
in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) 

generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")


# Use R packages for submap analysis  
#--------------------------------------#
# install.packages("C:\\Users\\NPC\\Documents\\WeChat Files\\wqwylyls\\FileStorage\\File\\2021-12\\GenePattern_1.0.2.tar.gz", type="source", repos=NULL)

library(GenePattern)
gctfileA <- 'submap_valid/GSE53819.for.SubMap.gct'
clsfileA <- 'submap_valid/GSE53819.for.SubMap.cls'

gctfileB <- 'submap_valid/GSE102349.for.SubMap.gct'
clsfileB <- 'submap_valid/GSE102349.for.SubMap.cls'


source('D:/GEO/script/submap.R')
submap.main(  gctfileA,
              gctfileB,
              clsfileA,
              clsfileB,
              output.filename="submap_valid/",
              ntag=100,
              nperm=50,
              nperm.fisher=1000,
              weighted.score.type=1,
              null.dist="pool",
              p.corr="Bonferroni",
              clust.row=1,
              clust.col=1,
              nom.p.mat="T",
              create.legend="T",
              rnd.seed=47365321) 
# Draw heatmap according to SubMap_SubMapResult.txt
library(pheatmap)
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#E7FFBF")
cherry    <- "#700353"
lightgrey <- "#dcddde"
mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

# Input the nominal p values in the file and the corrected P values in the heatmap
tmp <- matrix(c(0.987,0.412,0.997,0.001,0.040,0.388,0.086,0.976, 
                1,1,1,0.008,0.320,1,0.691,1), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("Cluster1_p","Cluster2_p","Cluster1_b","Cluster2_b"),c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")))

library(pheatmap)
submap.data.A <- matrix(c(0.002997003 , 0.981518482, 1.000000000 , 0.008991009),nrow = 2,byrow = T)
submap.data.B <- matrix(c(0.01198801 ,  1.00000000, 1.0000000, 0.03596404),nrow = 2,byrow = T)
submap.data <- rbind(submap.data.A,submap.data.B)
row.names(submap.data) <- c('Cluster1_p','Cluster2_p','Cluster1_b','Cluster2_b')
colnames(submap.data) <- c("GSE102349_Cluster1","GSE102349_Cluster2")
annotation_row <- data.frame(pvalue=rep(c('Nominal p value','Bonferroni corrected'),c(2,2)))
rownames(annotation_row) <- row.names(submap.data)
pdf(file = 'submap_valid/subclass_mapping.pdf',width= 10,height= 8)
pheatmap(submap.data,
         show_colnames = T,
         color = heatmap.YlGnPe[5:2],
         # display_numbers = matrix(ifelse(submap.data < 0.05,'p < .05',''),nrow(submap.data)),number_format = "%.3f",
         display_numbers = TRUE,
         number_format = "%.3f",
         cluster_rows = F,cluster_cols = F,
         annotation_row=annotation_row,
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         cellwidth=30,cellheight=30,main = "",fontsize_number = 9,
         #fontsize_row = 20,fontsize_col = 25,
         gaps_row=2,
         fontsize=12)
dev.off()
