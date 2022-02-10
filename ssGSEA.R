# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("GSVA")

library(GSVA)
library(ComplexHeatmap)
library(dplyr)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# read marker genes of immune cells
immunity <- read.csv("immune/immunitygene.csv", header = T)

# remove cell types that are not immune cells  
immunity <- immunity[!immunity$CellType %in% c("Blood vessels", "Normal mucosa", "SW480 cancer cells", "Lymph vessels"), ] %>% tidyr::separate_rows(Gene.Symbol,sep= " /// ")


immunity <- immunity %>% 
  split(., .$CellType) %>% 
  lapply(., function(x)(x$Symbol))
immunity <- lapply(immunity, unique)

# The infiltration level was quantified using ssGSEA

mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

expFile="01rawData/GSE102349_TPM_symbol.txt" 
expr=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)

expr[1:3,1:5]
class(expr)
expr1 <- as.matrix(expr)

gsva1 <- as.data.frame(t(gsva(expr1, immunity, method = "ssgsea")))

ssgseaOut <- gsva1
ssgseaOut_df=cbind(id=rownames(ssgseaOut),ssgseaOut)
ssgseaOut_df[1:3,1:6]
write.table(ssgseaOut_df,file="immune/GSE102349_ssGSEA_score.txt",sep="\t",quote=F,row.names=F)

gsva <- t(ssgseaOut)

library(pheatmap)
Cluster_res <- read.table("01rawData/GSE102349_Cluster_res.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

gsva <- gsva[,c(Cluster1,Cluster2)]
# Build comment information for the column
estimate <- read.table("immune/GSE102349_estimateScores.txt",header = T,sep = "\t",check.names = F,row.names = 1)
estimate <- estimate[c(Cluster1,Cluster2),]
annotation_col = data.frame(
  Type = factor(rep(c("Cluster1","Cluster2"), c(length(Cluster1) ,length(Cluster2))),levels = c("Cluster1","Cluster2")),
  StromalScore =estimate$StromalScore,
  ImmuneScore =estimate$ImmuneScore,
  ESTIMATEScore =estimate$ESTIMATEScore
)
rownames(annotation_col) = colnames(gsva)
head(annotation_col)
# A list of colors for custom comment information
ann_colors = list(
  Type = c(Cluster1 = mycol[2], Cluster2 = mycol[3] )
 )
head(ann_colors)


gsva <- gsva[apply(gsva, 1, function(x) sd(x)!=0),]

pdf(file="immune/GSE102349_ssGSEA_score_heatmap.pdf",height=8,width=10)
pheatmap(gsva,
         scale = "row", 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F, 
         cluster_cols =F,
         border=FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dev.off()


# draw a boxplot
boxdata <-  as.data.frame(t(gsva))
risk <- read.table("01rawData/GSE102349_Cluster_res.txt",header = T,sep = "\t",check.names = F)

boxdata$id <- rownames(boxdata)
boxdata <- merge(boxdata,risk[,c("id","Cluster")],by="id")
write.table(boxdata, "immune/GSE102349_ssGSEA_score_cluster.txt", quote = F, row.names = F,sep = "\t")

boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)
boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])

p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
               color = "Cluster", palette = mycol[c(2:3,8)],
               add = "jitter",
               add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("ssGSEA score")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")
# p2 + stat_compare_means(comparisons = split(t(combn(levels(boxdata1$Cluster),2)),1:nrow(t(combn(levels(boxdata1$Cluster),2)))))
ggsave("immune/GSE102349_ssGSEA_score_boxplot.pdf",height=5,width=10)


tcga_expr <- expr
tcga_gsva <-as.data.frame(t(gsva))

genelist <- c("RAP1GAP","GNL2","PSD4","COBLL1","CDV3","NR1H3","G6PC")
immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="pearson")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}


# Using FOXP3 as an example, test the function
immuscore("RAP1GAP")
# The correlation between genelist and immune infiltration was calculated in batches
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
# Save to file
write.csv(data, "06ssGSEA/immune_cells_modelGenes_cor.csv", quote = F, row.names = F)

data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]

ggplot(data, aes( immune_cells,gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),# Don't the title  
        axis.ticks.x=element_blank(),# Don't x
        axis.title.y=element_blank(),# Don't y
        axis.text.x = element_text(angle = 45, hjust = 1),# Adjust the X-axis text
        axis.text.y = element_text(size = 8))+# Adjust the Y-axis text
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
ggsave("06ssGSEA/immune_cells_modelGenes_correlation.pdf", width = 8.5, height = 4)
