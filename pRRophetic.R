
rm(list = ls())  
options(stringsAsFactors = F)
expFile="01rawData/GSE102349_TPM_symbol.txt"                                                 
pdFile = "01rawData/GSE102349_Cluster_res.txt"

testExpr <- read.table(expFile,sep="\t",header=T,row.names=1,check.names=F)
pd<- read.table(pdFile,sep="\t",header=T,row.names = 1,check.names=F)

# drug prediction R package of pRRophetic
library(pRRophetic)
library(ggplot2)
library(ggpubr)
library(cowplot)
mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

# drug name
GCP.drug <- read.table("pRRophetic/drug_sig.txt") 
GCP.drug <- GCP.drug$V1

exprData <- as.matrix(testExpr[,rownames(pd)])
# customize enough box colors that the number of colors is at least 
# equal to the number of groups
jco <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
         '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

# drug sensitivity prediction
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # Initialization list
plotp <- list()

for (drug in GCP.drug) {
  set.seed(1248103) # because the prediction process defaults to a 10-fold CV, the seeds are set so that the results are repeatable
  cat(drug," starts!\n") # indicates that the current drug has been analyzed
  
  # To predict the IC50 value, use the default parameters. 
  # For details, see pRRopheticPredict parameter list
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = exprData ,
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              batchCorrect = "eb",
                                              selection = 1, 
                                              dataset = "cgp2014")
  if(!all(names(predictedPtype[[drug]])==rownames(pd))) {stop("Name mismatched!\n")} 
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "Cluster"=ifelse(pd$Cluster == "Cluster1","Cluster1","Cluster2"), 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$Cluster <- factor(predictedBoxdat[[drug]]$Cluster,levels = c("Cluster1","Cluster2"),ordered = T) 
  # drawing
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=Cluster, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = Cluster),outlier.colour = NA,notch = T,size = 0.3)+
    geom_jitter(aes(fill = Cluster),shape = 21,size=2,width = 0.2)+
    geom_violin(aes(fill = Cluster),position = position_dodge(width = .75), 
                size = 0.3,alpha = 0.4,trim = T)+
    scale_fill_manual(values=mycol[c(2,3)]) + 
    theme_classic()+
    stat_compare_means(comparisons = list(c('Cluster1','Cluster2')),method="wilcox.test", label = "p.format")+
    theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) # add the title

  plotp[[drug]] <- p 
  cat(drug," has been finished!\n") 
}

# Merged images
# Show the first four sensitive drugs
p1 <- plot_grid(plotp[[1]],plotp[[2]],plotp[[3]],plotp[[4]],labels = c("A","B","C","D"),nrow = 1) # title可以AI下拉到合适位置，就如例文所示
ggsave("pRRophetic/TOP4 boxplot of predicted IC50.pdf", width = 8, height = 3.6)

# display multiple sensitive drugs
p2 <- plot_grid(plotlist=plotp, ncol=5)
ggsave("pRRophetic/boxplot of predicted IC50_multiple.pdf", width = 10, height = 15)

# The differences between groups were examined
p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$Cluster %in% "Cluster1"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$Cluster %in% "Cluster2"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # The rank sum test for p values of the two groups of samples
}
names(p) <- GCP.drug
print(p) 
# Save to file
write.table(p,"pRRophetic/output_pvalue.txt", quote = F, sep = "\t")

# TIDE was used to compare the immune responses of different subtypes
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
TIDE <- testExpr

TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"pRRophetic/TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)


TIDE.res <- read.csv("pRRophetic/TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
pd$TIDE <- TIDE.res[rownames(pd),"Responder"]
print(table(pd$TIDE,pd$Cluster))
# To test whether immunotherapy response is related to subtypes, P<0.05 is relevant
print(fisher.test(table(pd$TIDE,pd$Cluster))) 

# submap predict immunotherapy response of subtypes
# custom functions are used to generate the data format required by submap
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

# The data format needed to create the submap (SKCM)
skcm.immunotherapy.logNC <- read.table("pRRophetic/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) 
skcm.immunotherapy.info <- read.table("pRRophetic/skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

# The data format needed to create the submap  (GEO)
tmp <- read.table(expFile,sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1) 
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # take the list of genes after intersection

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# The file name that generates the output data
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")


samples.C1 <- rownames(pd[which(pd$Cluster == "Cluster1"),])
samples.C2 <- rownames(pd[which(pd$Cluster == "Cluster2"),])

sam_info <- data.frame("Cluster"=c(samples.C1,samples.C2),row.names = c(samples.C1,samples.C2))
sam_info$rank <- rep(c(1,2),times=c(length(samples.C1),length(samples.C2))) 

# The file name that generates the output data
gct_file <- "GSE102349.Immune2.for.SubMap.gct"
cls_file <- "GSE102349.Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) 
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# GenePattern was used for submap analysis

library(GenePattern)
gctfileA <- 'pRRophetic/GSE102349.Immune2.for.SubMap.gct'
clsfileA <- 'pRRophetic/GSE102349.Immune2.for.SubMap.cls'

gctfileB <- 'pRRophetic/skcm.immunotherapy.for.SubMap.gct'
clsfileB <- 'pRRophetic/skcm.immunotherapy.for.SubMap.cls'

source('D:/GEO/script/submap.R')
submap.main(  gctfileA,
              gctfileB,
              clsfileA,
              clsfileB,
              output.filename="pRRophetic/SubMap",
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

library(pheatmap)
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"


# Input the nominal p values in the file and the corrected P values in the heat map
tmp <- matrix(c(0.987,0.412,0.997,0.001,0.040,0.388,0.086,0.976, 
                1,1,1,0.008,0.320,1,0.691,1), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("Cluster1_p","Cluster2_p","Cluster1_b","Cluster2_b"),c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")))

library(pheatmap)
submap.data.A <- matrix(c(0.9700300 , 0.4595405 , 0.9930070 , 0.000999001,
                          0.1878122 , 0.3456543 , 0.1678322 , 0.955044955),nrow = 2,byrow = T)
submap.data.B <- matrix(c(1,  1, 1, 0.007992008,
                          1,  1, 1, 1.000000000),nrow = 2,byrow = T)
submap.data <- rbind(submap.data.A,submap.data.B)
row.names(submap.data) <- c('Cluster1_p','Cluster2_p','Cluster1_b','Cluster2_b')
colnames(submap.data) <- c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")
annotation_row <- data.frame(pvalue=rep(c('Nominal p value','Bonferroni corrected'),c(2,2)))
rownames(annotation_row) <- row.names(submap.data)
pdf(file = 'pRRophetic/subclass mapping.pdf',width= 12,height= 8)
pheatmap(submap.data,
         show_colnames = T,
         color = heatmap.YlGnPe[5:1],
         display_numbers = matrix(ifelse(submap.data < 0.05,'p < .05',''),nrow(submap.data)),number_format = "%.3f",
         cluster_rows = F,cluster_cols = F,
         annotation_row=annotation_row,
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         cellwidth=30,cellheight=30,main = "",fontsize_number = 9,
         #fontsize_row = 20,fontsize_col = 25,
         gaps_row=2,
         fontsize=12)
dev.off()


# Draw a boxplot of the TIDE response
ann <-read.table("01rawData/GSE102349_Cluster_res.txt",header = T,sep = "\t",check.names = F)
head(ann)
print(table(ann$Cluster))

TIDE.res <- read.csv("pRRophetic/TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)

TIDE.res$id <- rownames(TIDE.res)

# To test whether immunotherapy response is related to subtypes, P<0.05 is relevant
ann$TIDE <- TIDE.res[ann$id,"Responder"]
print(table(ann$TIDE,ann$Cluster))
print(fisher.test(table(ann$TIDE,ann$Cluster))) 



boxdata <- merge(TIDE.res,ann,by="id")
boxdata$Cluster <- factor(boxdata$Cluster,levels=c("Cluster1","Cluster2"))


ggplot(boxdata,aes(Cluster,TIDE.x,fill=Cluster))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+ ylab('TIDE')+ggtitle('TIDE prediction')+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(2,3)])+
  stat_compare_means(comparisons = list(c('Cluster1','Cluster2')),
                     label = 'p.signif')+
  stat_compare_means(label.y = max(boxdata$TIDE.x)+0)
ggsave("pRRophetic/TIDE_boxplot.pdf",width = 3.5,height = 3.5)
