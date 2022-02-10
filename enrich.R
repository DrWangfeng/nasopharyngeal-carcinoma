rm(list = ls())
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("GOSemSim")
library("enrichplot")
library("ggplot2")
                 
rt=read.table("DEGs/GSE102349_edgeR_diff.txt",sep="\t",header=T,check.names=F,row.names = 1)

entrezIDs <- mget(rownames(rt), org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out = cbind(symbol=rownames(rt),entrezID=entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
out = out[,c(1,6)]
out = cbind(symbol=rownames(out),out)
write.table(out,file="DEGs/id.txt",quote=F,row.names = F,sep="\t")

rt <- read.table("DEGs/id.txt",sep="\t",header = T)
entrezID_gene <- na.omit(rt$entrezID)

x <- enrichGO(gene = entrezID_gene,
               OrgDb = org.Hs.eg.db,
               pvalueCutoff =0.05,
               qvalueCutoff = 0.05,
               ont="all",
               readable =T
)
write.table(x,file="DEGs/GO.txt",quote=F,row.names = F,sep = "\t")


pdf(file="DEGs/GO_barplot.pdf",width = 9.5,height = 9)
barplot(x, drop = TRUE, showCategory =10,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
dev.off()


pdf(file="DEGs/GO_bubble.pdf",width = 9.5,height = 9)
dotplot(x,showCategory = 10,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
dev.off()

kk <- enrichKEGG(gene = entrezID_gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)

y <- setReadable(kk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.table(y,file="DEGs/KEGG.txt",quote=F,row.names = F,sep = "\t")

pdf(file="DEGs/KEGG_barplot.pdf",width = 9,height = 4)
barplot(kk, drop = TRUE, showCategory = 10)
dev.off()

pdf(file="DEGs/KEGG_bubble.pdf",width = 9,height = 4)
dotplot(kk, showCategory = 10)
dev.off()
