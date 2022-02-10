# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("GSVA")

library(GSVA)
library(GSEABase)
library(limma)
library(tidyverse)
library(clusterProfiler)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

exprFile <- "01rawData/GSE102349_TPM_symbol.txt" 
gmtFile <- "GSVA/c2.cp.kegg.v7.4.symbols.gmt"
expr <- read.table(exprFile,sep="\t",header=T,check.names=F,row.names = 1)
# # ~~~~~~~~~~~
# gmt <- read.table(gmtFile,sep="\t",check.names=F,row.names = 1)
# geneSetlist <- list()
# geneSet <- as.data.frame(t(gmt))
# for (i in colnames(geneSet)) {
#   geneSetlist[[i]] <-  geneSet[,i][which(geneSet[,i]!="")]
# 
# }
# geneSet <- geneSetlist

geneSet=getGmt(gmtFile, 
                collectionType=BroadCollection(category="c2"),
               geneIdType=SymbolIdentifier())
# The default value is "Gaussian".  The parrallel. Sz parameter sets the number of 
# parallel threads, because each gene set is evaluated independently.  

gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
               min.sz=5, 
               max.sz=500, 
               kcdf="Gaussian", 
               method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="GSVA/GSVA_KEGG_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)


# kruskal.test Test to test pathway differences
samples_input <- "01rawData/GSE102349_Cluster_res.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)

gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file, "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  # kruskalTest <- kruskal.test(get(pathway)~Cluster,gsva_result_df)
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster2Med=median(rt[rt$Cluster=="Cluster2","expression"])
  Cluster1Med=median(rt[rt$Cluster=="Cluster1","expression"])
  diffMed=Cluster1Med-Cluster2Med
  
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1Med,
                        Cluster2Mean=Cluster2Med,
                        logFC=diffMed,
                        p.Value=pvalue))
  }
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "GSVA/GSVA_KEGG_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)



gmtFile <- "GSVA/h.all.v7.4.symbols.gmt"
geneSet=getGmt(gmtFile, 
               collectionType=BroadCollection(category="h"),
               geneIdType=SymbolIdentifier())

gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
                    min.sz=5, 
                    max.sz=500, 
                    kcdf="Gaussian", 
                    method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="GSVA/GSVA_HALLMARK_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)


# kruskal.test Test to test pathway differences
samples_input <- "01rawData/GSE102349_Cluster_res.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)
gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file, "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster2Med=median(rt[rt$Cluster=="Cluster2","expression"])
  Cluster1Med=median(rt[rt$Cluster=="Cluster1","expression"])
  diffMed=Cluster1Med-Cluster2Med
  
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1Med,
                        Cluster2Mean=Cluster2Med,
                        logFC=diffMed,
                        p.Value=pvalue))
}
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "GSVA/GSVA_HALLMARK_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)


gmtFile <- "GSVA/c5.go.bp.v7.4.symbols.gmt"
geneSet=getGmt(gmtFile, 
               collectionType=BroadCollection(category="c5"),
               geneIdType=SymbolIdentifier())

gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
                    min.sz=5, 
                    max.sz=500, 
                    kcdf="Gaussian", 
                    method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="GSVA/GSVA_GOBP_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)

# kruskal.test test to test pathway differences
samples_input <- "01rawData/GSE102349_Cluster_res.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)

gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file, "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster2Med=median(rt[rt$Cluster=="Cluster2","expression"])
  Cluster1Med=median(rt[rt$Cluster=="Cluster1","expression"])
  diffMed=Cluster1Med-Cluster2Med
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1Med,
                        Cluster2Mean=Cluster2Med,
                        logFC=diffMed,
                        p.Value=pvalue))
}
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "GSVA/GSVA_GOBP_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)


gmtFile <- "GSVA/c2.cp.reactome.v7.4.symbols.gmt"
geneSet=getGmt(gmtFile, 
               collectionType=BroadCollection(category="c2"),
               geneIdType=SymbolIdentifier())

gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
                    min.sz=5, 
                    max.sz=500, 
                    kcdf="Gaussian", 
                    method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="GSVA/GSVA_reactome_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)


samples_input <- "01rawData/GSE102349_Cluster_res.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)

gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file, "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster2Med=median(rt[rt$Cluster=="Cluster2","expression"])
  Cluster1Med=median(rt[rt$Cluster=="Cluster1","expression"])
  diffMed=Cluster1Med-Cluster2Med
  
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1Med,
                        Cluster2Mean=Cluster2Med,
                        logFC=diffMed,
                        p.Value=pvalue))
}
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "GSVA/GSVA_reactome_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# GSVA enrichment analysis was performed on the gene set obtained by ourselves
gmtFile <- "GSVA/GSVA_geneSet.txt"

gmt <- read.table(gmtFile,sep="\t",check.names=F,row.names = 1)
geneSetlist <- list()
geneSet <- as.data.frame(t(gmt))
for (i in colnames(geneSet)) {
  geneSetlist[[i]] <-  geneSet[,i][which(geneSet[,i]!="")]

}
geneSet <- geneSetlist

# The default value is "Gaussian".  The parrallel. sz parameter sets 
# the number of parallel threads, because each gene set is evaluated independently.  
gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
                    min.sz=1, 
                    max.sz=500, 
                    kcdf="Gaussian", 
                    method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="GSVA/GSVA_geneSet_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)

# Kruskal.test to test pathway differences
samples_input <- "01rawData/GSE102349_Cluster_res.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)

gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file, "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  # kruskalTest <- kruskal.test(get(pathway)~Cluster,gsva_result_df)
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster2Med=median(rt[rt$Cluster=="Cluster2","expression"])
  Cluster1Med=median(rt[rt$Cluster=="Cluster1","expression"])
  diffMed=Cluster1Med-Cluster2Med
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1Med,
                        Cluster2Mean=Cluster2Med,
                        logFC=diffMed,
                        p.Value=pvalue))
}
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "GSVA/GSVA_geneSet_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# compare custom pathway differences
# Draw customized gene sets for GSVA heatmap of different m6A subtypes
exprSet_m6A <-  read.table("GSVA/GSVA_geneSet_result.txt",header = T,sep = "\t",check.names = F,row.names = 1)
outDiff.mRNA <-  read.table("GSVA/GSVA_geneSet_wilcoxTest.txt",header = T,sep = "\t",check.names = F,row.names = 1)
geneNum=30
# outDiff.mRNA=outDiff.mRNA[outDiff.mRNA$p.adjust < 0.05,]
outDiff.mRNA=outDiff.mRNA[order(as.numeric(as.vector(outDiff.mRNA$logFC))),]
diffGeneName=rownames( outDiff.mRNA )
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet_m6A[hmGene,]

max(hmExp)
min(hmExp)

library(pheatmap)
Cluster_res <- read.table("01rawData/GSE102349_Cluster_res.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]


hmExp <-as.matrix(hmExp[,c(Cluster1,Cluster2)])

# build column comment information
annotation_col = data.frame(
  Type = factor(rep(c("Cluster1","Cluster2"), c(length(Cluster1) ,length(Cluster2))),levels = c("Cluster1","Cluster2"))
)
rownames(annotation_col) = colnames(hmExp)
head(annotation_col)
# 自定注释信息的颜色列表 "#E69F00","#56B4E9"
ann_colors = list(
  Type = c(Cluster1 = mycol[2], Cluster2 = mycol[3] )
)
head(ann_colors)


loc <- order(annotation_col$Type,colSums(hmExp),decreasing = F)


pdf(file="GSVA/GSVA_geneSet_heatmap.pdf",height=5,width=9)
pheatmap(#hmExp,
  hmExp[,loc],
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = F, 
  cluster_cols =F,
  border=FALSE,
  cellwidth = 4,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
scores =read.table("GSVA/GSVA_geneSet_result.txt",header=T,sep="\t",check.names=F,row.names = 1)
boxdata <- as.data.frame(t(scores)) %>% rownames_to_column("id")
risk <- read.table("01rawData/GSE102349_Cluster_res.txt",header = T,sep = "\t",check.names = F)

boxdata <- merge(boxdata,risk[,c("id","Cluster")],by="id")
write.table(boxdata, "GSVA/GSVA_geneSet_result_cluster.txt", quote = F, row.names = F,sep = "\t")


boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)
mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
# boxdata1$risk <- factor(boxdata1$risk,levels = c("Low Risk","High Risk"))
# boxdata1$Cluster <- factor(boxdata1$Cluster,levels = c("Cluster1","Cluster2"))
p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
                color = "Cluster", palette = mycol[c(2:3,8)],
                add = "jitter",
                add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("ssGSEA score")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")
# p2 + stat_compare_means(comparisons = split(t(combn(levels(boxdata1$Cluster),2)),1:nrow(t(combn(levels(boxdata1$Cluster),2)))))
ggsave("GSVA/GSVA_geneSet_ssGSEA_score_boxplot.pdf",height=5.5,width=11)


# the differential genes of two cluters
library(edgeR)
exprFile <- "01rawData/GSE102349_counts_symbol.txt"
expr <- read.table(exprFile,sep="\t",header=T,check.names=F,row.names = 1)
expr[1:4,1:7]
samples_input <- "01rawData/GSE102349_Cluster_res.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)
samples_file <- samples_file[order(samples_file$Cluster),]
samples_file[1:5,1:2]

Cluster1 <- samples_file$id[samples_file$Cluster =="Cluster1"]
Cluster2 <- samples_file$id[samples_file$Cluster =="Cluster2"]
expr <- expr[,c(Cluster1,Cluster2)]

# Construct a DGEList object
dgelist <- DGEList(counts = expr, group = samples_file$Cluster)

#
keep <- rowSums(cpm(dgelist) > 1 ) >= 10
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dim(dgelist)
# Filtering low count data, such as CPM standardization
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')


design <- model.matrix(~0+ samples_file$Cluster)
rownames(design) <- samples_file$id
colnames(design) <- c("Cluster1" ,"Cluster2")

# Estimate the dispersion of gene expression values
dge <- estimateDisp(dgelist_norm, design)
# Negative binomial generalized log-linear model
fit <- glmFit(dge, design)
lrt <- glmLRT(fit,contrast=c(1,-1)) 
result <- topTags(lrt, n = nrow(dgelist$counts))
allDiff <- result$table

logFCCutoff <- 1
pvalueCutoff <- 0.05

outDiff=allDiff[(abs(allDiff$logFC)>logFCCutoff & allDiff$FDR<pvalueCutoff),]
outDiff <- outDiff %>% rownames_to_column(var = "id")
write.table(outDiff,file="DEGs/GSE102349_edgeR_diff.txt",row.names=F,quote=F,sep = "\t")

allDiff <- allDiff  %>% rownames_to_column(var = "id")
write.table(allDiff, 'DEGs/GSE102349_edgeR_alldiff.txt', sep = '\t', row.names=F, quote = FALSE)
