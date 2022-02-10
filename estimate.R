# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")
rm(list=ls())

library(limma)
library(estimate)

inputFile="01rawData/GSE102349_TPM_symbol.txt"

# estimate
filterCommonGenes(input.f= inputFile,
                  output.f="immune/GSE102349_commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "immune/GSE102349_commonGenes.gct",
              output.ds="immune/GSE102349_estimateScore.gct", 
              platform="illumina")


scores=read.table("immune/GSE102349_estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(id=colnames(scores),scores)
write.table(out,file="immune/GSE102349_estimateScores.txt",sep="\t",quote=F,col.names=F)



# Draw the boxplot of estimate score in the cluster grouping
library(tidyr)
library(ggpubr)

mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

risk <- read.table("01rawData/GSE102349_Cluster_res.txt",header = T,sep = "\t",check.names = F)
scores =read.table("immune/GSE102349_estimateScores.txt",header=T,sep="\t",check.names=F)
Data <- merge(risk,scores,by="id")


plotData <- Data[,c("Cluster","StromalScore","ImmuneScore","ESTIMATEScore")]
plotData$Cluster <- factor(plotData$Cluster,levels = c("Cluster1","Cluster2"))

p1 <- ggplot(plotData,aes(Cluster,StromalScore,fill=Cluster))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(2,3)])+
  stat_compare_means(comparisons = list(c('Cluster1','Cluster2')),
                      label = 'p.signif')+
  stat_compare_means(label.y = max(plotData$StromalScore)+5.5)

p2 <-  ggplot(plotData,aes(Cluster,ImmuneScore ,fill=Cluster))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(2,3)])+
  stat_compare_means(comparisons = list(c('Cluster1','Cluster2')),
                     label = 'p.signif')+
  stat_compare_means(label.y = max(plotData$ImmuneScore )+5.5)

p3 <- ggplot(plotData,aes(Cluster,ESTIMATEScore ,fill=Cluster))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(2,3)])+
  stat_compare_means(comparisons = list(c('Cluster1','Cluster2')),
                     label = 'p.signif')+
  stat_compare_means(label.y = max(plotData$ESTIMATEScore )+5.5)
# Group presentation on one page using the ggpubr package's function ggarrange ()  
ggarrange(p1,p2,p3,labels = c("A", "B","C"),ncol = 3, nrow = 1)

ggsave("immune/GSE102349_Score_violin.pdf",width = 11,height = 3.5)

