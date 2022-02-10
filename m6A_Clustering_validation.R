rm(list=ls())
# Load the R package
# BiocManager::install("paletteer")
library(ConsensusClusterPlus)
library(tidyverse)
library(pheatmap)
library(survival)
library(survminer)
library(Rtsne)
library(ggplot2)
# Define drawing colors
mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

# Read the expression data of GSE102349 and do univariate COX analysis of m6A gene
exprSet0 <- read.table("01rawData/GSE102349_TPM_symbol.txt",row.names = 1,check.names = F,sep="\t",header = T)
HYPOXIA <- read.table("01rawData/m6A_genes.txt",header = T,sep="\t")

clincal <- read.table("01rawData/GSE102349_clinical_PFS.txt",sep = "\t",header = T,check.names = F)

exprSet_m6A <- exprSet0[rownames(exprSet0)  %in% HYPOXIA$gene,clincal$id]
exprSet_m6A1 <- as.data.frame( t(exprSet_m6A))
# exprSet$Id <- rownames(exprSet)
exprSet_m6A1$id <- substr(rownames(exprSet_m6A1),1,12)
exprSet_m6A1_Time <- merge(clincal,exprSet_m6A1,by="id")
exprSet_m6A1_Time[1:3,1:4]

write.table(exprSet_m6A1_Time,file = "03Cox/GSE102349_m6A_exprSetTime.txt",row.names = F,quote = F,sep = "\t")
# univariate COX analysis of m6A gene
library(survival)
rt=read.table("03Cox/GSE102349_m6A_exprSetTime.txt",header=T,check.names=F,row.names=1,sep = "\t") 
rt1 <- rt
rt$event <- ifelse(rt$event =="Last follow-up",0,1)
sur<-Surv(time=rt$`time to event`, event = rt$event)

pFilter=1
outTab=data.frame()
sigGenes = colnames(rt1)[1:10]
# Select the name of the variable to be used for univariate analysis
(variable_names<-colnames(rt)[-c(1:10)])

for(gene in variable_names){
  if(sd(rt[,gene])<0.01){next}
  if(grepl("-", gene)){next}
  # cox=coxph(as.formula(paste0('sur~',gene))  ,data = rt)
  cox=coxph(Surv(`time to event`, event) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
    diff=survdiff(sur ~group,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    if(pValue<pFilter){
      sigGenes=c(sigGenes,gene)
      outTab=rbind(outTab,
                   cbind(gene=gene,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         coxPvalue=coxP) )
    }
  }
}

write.table(outTab,file="03Cox/GSE102349_m6A_uniCox.txt",sep="\t",row.names=F,quote=F)   


HYPOXIA <- read.table("03Cox/GSE102349_m6A_uniCox.txt",header = T,sep="\t")
m6A_gene <- HYPOXIA$gene[HYPOXIA$coxPvalue < 0.05]
exprSet_m6A <- exprSet0[rownames(exprSet0)  %in% m6A_gene,clincal$id]
feExpr_matrix <- as.matrix(log2(exprSet_m6A+1) )

pheatmap(feExpr_matrix)

# The results will output the classification of k from 2-6 in each case. 
# The clustering method is HC with a sampling ratio of 0.8. PNG images 
# will be output at last  
Cluster  <- ConsensusClusterPlus(feExpr_matrix, 
                                 maxK = 6, # Maximum number of clusters
                                 reps = 1000, # Number of resamples 
                                 pItem = 0.8, # Percentage of samples to be re-sampled
                                 pFeature = 1,  
                                 clusterAlg = "hc", # The clustering algorithm used
                                 distance = "euclidean",  # method of calculating distance
                                 seed = 99,
                                 title = "02ConsensusCluster",
                                 plot = "png")


maxK = 6
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = Cluster[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i
# The optimal K
optK = Kvec[which.min(PAC)]
optK

# Extract the clustering results
annCol <- data.frame(Cluster = paste0("Cluster",
                                      Cluster[[optK]][["consensusClass"]]),
                     row.names = colnames(feExpr_matrix))

Cluster_res <- annCol %>%
  as.data.frame() %>%
  rownames_to_column("id")# %>% arrange(Cluster)
table(Cluster_res$Cluster)
write.table(Cluster_res,file = "01rawData/GSE102349_Cluster_res.txt",row.names = F,quote = F,sep = "\t")

Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

exprSet_m6A <- exprSet_m6A[,c(Cluster1,Cluster2)]
exprSet_m6A_out <- cbind(gene_name=rownames(exprSet_m6A),exprSet_m6A)
exprSet_m6A_out[1:3,1:5]
write.table(exprSet_m6A_out,file = "03Cox/GSE102349_m6A_exprSet.txt",row.names = F,quote = F,sep = "\t")



# Loading tsne package  
library(Rtsne)
set.seed(99)
# T-SNE dimension reduction analysis was performed using Rtsne function  
feExpr1 <- as.matrix(t(feExpr_matrix))
tsne_out <- Rtsne(feExpr1,pca=FALSE,dims=2,
                  perplexity=15,theta=0.0) # Run TSNE
plot(tsne_out$Y,col=factor( Cluster_res$Cluster),asp=1)
tsne_plot  = data.frame(tSNE1  =tsne_out$Y[,1], tSNE2  = tsne_out$Y[,2],Cluster = Cluster_res$Cluster)

head(tsne_plot)
# The tSNE dimension reduction cluster diagram was drawn
ggplot(tsne_plot,aes(x=tSNE1,y=tSNE2))+
  geom_point(aes(fill = Cluster),shape = 21)+
  stat_ellipse(aes(color = Cluster,fill = Cluster),
               geom = "polygon",
               level = 0.95,
               alpha = 0.3,
               linetype = 2)+
  scale_color_manual(values = mycol[c(2,3,8)])+
  scale_fill_manual(values= mycol[c(2,3,8)])+
  theme_classic()+
  theme(legend.position = "top")

ggsave("02ConsensusCluster/GSE102349_Cluster_tSNE.pdf",width = 4.2,height = 4.5)


Cluster_res <- Cluster_res %>%
  inner_join(clincal, "id")
head(Cluster_res)
Cluster_res$event <- ifelse(Cluster_res$event =="Last follow-up",0,1)
# KM survival analysis 
fit <- survfit(Surv(`time to event`, `event`) ~ Cluster, 
               data = Cluster_res)
ggsurvplot(fit,
           pval = TRUE,
           linetype = "solid",  
           palette = mycol[c(2,3,8)],
           #surv.median.line = "hv", 
           title = "GSE102349",
           ylab = "Progression-free survival (percentage)",
           xlab = " Time (Years)",
           # legend.title = "Survival Plot",
           legend = c(0.8,0.30),
           # legend.labs = c("Cluster1","Cluster2","Cluster3"),
           legend.labs = c("Cluster1","Cluster2"),
           legend.title="",
           risk.table = T,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)
dev.copy2pdf(file = "02ConsensusCluster/GSE102349_Cluster_PFS_survival.pdf", width = 5,height = 5.5)
dev.off()


# The external GSE53819 dataset validates the tSNE clustering
data <- read_tsv("01rawData/GSE53819_NPC_exprSet.txt",)
data_valid <-  data[data$symbol %in% m6A_gene, ] %>% column_to_rownames("symbol")
feExpr_matrix <- as.matrix((data_valid) )

Cluster  <- ConsensusClusterPlus(feExpr_matrix, 
                                 maxK = 6, # Maximum number of clusters
                                 reps = 1000, # Number of resamples 
                                 pItem = 0.8, # Percentage of samples to be re-sampled
                                 pFeature = 1,  
                                 clusterAlg = "hc", # The clustering algorithm used
                                 distance = "euclidean",  # method of calculating distance
                                 seed = 99,
                                 title = "02ConsensusCluster_valid",
                                 plot = "png")


maxK = 6
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = Cluster[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i
# The optimal K
optK = Kvec[which.min(PAC)]
optK

# optK=3
# # Extract the clustering results
annCol <- data.frame(Cluster = paste0("Cluster",
                                      Cluster[[optK]][["consensusClass"]]),
                     row.names = colnames(feExpr_matrix))
head(annCol)


Cluster_res <- annCol %>%
  as.data.frame() %>%
  rownames_to_column("id")# %>% arrange(Cluster)
table(Cluster_res$Cluster)
write.table(Cluster_res,file = "01rawData/GSE53819_Cluster_res.txt",row.names = F,quote = F,sep = "\t")

Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

data_valid_m6A <- data_valid[,c(Cluster1,Cluster2)]
data_valid_m6A_out <- cbind(gene_name=rownames(data_valid_m6A),data_valid_m6A)
data_valid_m6A_out[1:3,1:5]
write.table(data_valid_m6A_out,file = "03Cox/GSE53819_NPC_m6A_exprSet.txt",row.names = F,quote = F,sep = "\t")




library(Rtsne)
set.seed(99)
# T-SNE dimension reduction analysis was performed using Rtsne function 
feExpr1 <- as.matrix(t(feExpr_matrix))
tsne_out <- Rtsne(feExpr1,pca=FALSE,dims=2,
                  perplexity=3,theta=0.0) # Run TSNE
plot(tsne_out$Y,col=factor( Cluster_res$Cluster),asp=1)
tsne_plot  = data.frame(tSNE1  =tsne_out$Y[,1], tSNE2  = tsne_out$Y[,2],Cluster = Cluster_res$Cluster)

head(tsne_plot)

mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")
# The tSNE dimension reduction cluster diagram was drawn
ggplot(tsne_plot,aes(x=tSNE1,y=tSNE2))+
  geom_point(aes(fill = Cluster),shape = 21)+
  stat_ellipse(aes(color = Cluster,fill = Cluster),
               geom = "polygon",
               level = 0.9,
               alpha = 0.3,
               linetype = 2)+
  scale_color_manual(values = mycol[c(2,3,8)])+
  scale_fill_manual(values= mycol[c(2,3,8)])+
  theme_classic()+
  theme(legend.position = "top")

ggsave("02ConsensusCluster_valid/GSE53819_Cluster_tSNE.pdf",width = 4.2,height = 4.5)

