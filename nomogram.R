library(survival)
library(regplot)
library(rms)
library(nomogramEx)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
pbc=read.table("model/wholeSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F)       #??ȡ?????ļ?

pbc$event <- ifelse(pbc$event =="Last follow-up",0,1)
(colnames(pbc)[10:ncol(pbc)-2])
dd <- datadist(pbc)
options(datadist="dd")
options(na.action="na.delete")

units(pbc$`time to event`) <- "Month" 
summary(pbc$`time to event`)
coxpbc <- cph(formula = Surv(`time to event`,event) ~  REEP2 +TMSB15A+DSEL+ID4 ,data=pbc,x=T,y=T,surv = T,na.action=na.delete)  #,time.inc =2920

print(coxpbc)

surv <- Survival(coxpbc) 
surv1 <- function(x) surv(12,x)
surv2 <- function(x) surv(24,x)
surv3 <- function(x) surv(36,x)


x <- nomogram(coxpbc,fun = list(surv1,surv2,surv3),lp=T,
              funlabel = c('1-year survival Probability','2-year survival Probability','3-year survival Probability'),
              maxscale = 100,fun.at = c(0.9,0.85,0.8,0.7,0.5,0.3,0.1))

pdf("nomogram/nomogram_classical.pdf",width = 11,height = 7)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.2,cex.axis = 0.8,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

# Calibration Curve is drawn for verification
f1 <- cph(formula =  Surv(`time to event`,event) ~REEP2 +TMSB15A+DSEL+ID4 ,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 12) 
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=12,m=18,B=500) 

f3 <- cph(formula =  Surv(`time to event`,event) ~ REEP2 +TMSB15A+DSEL+ID4 ,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 24) 
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=24,m=18,B=500) 

f5 <- cph(formula =  Surv(`time to event`,event) ~ REEP2 +TMSB15A+DSEL+ID4 ,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 36) 
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=36,m=18,B=500) 

pdf("nomogram/calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", 
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced PFS (%)",ylab = "Observed PFS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#F0027F"),
     xlim = c(0,1),ylim= c(0,1),col = c("#F0027F"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#F0027F"), pch = 16)

mtext("")

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)


abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", 
       legend = c("1-year","2-year","3-year"),
       col =c("#2166AC","#F0027F","#B2182B"), 
       lwd = 2,
       cex = 1.2,
       bty = "n")
dev.off()



pdf("nomogram/calibration_1y.pdf",width = 8,height = 8)
plot(cal1,
     lwd = 2,
     lty = 1,
     errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced PFS (%)",ylab = "Observed PFS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) 
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', 
      lwd = 2, 
      pch = 16, 
      col = c("#2166AC")) 
mtext("")
box(lwd = 1) 
abline(0,1,lty = 3, 
       lwd = 2, 
       col = c("#224444")
) 
dev.off()

pdf("nomogram/calibration_2y.pdf",width = 8,height = 8)
plot(cal3,
     lwd = 2,
     lty = 1,
     errbar.col = c("#F0027F"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced PFS (%)",ylab = "Observed PFS (%)",
     col = c("#F0027F"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],
      type= 'b',
      lwd = 2,
      col = c("#F0027F"),
      pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,
       lwd = 2,
       col =c("#224444"))
dev.off()


pdf("nomogram/calibration_3y.pdf",width = 8,height = 8)
plot(cal5,
     lwd = 2,
     lty = 1,
     errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced PFS (%)",ylab = "Observed PFS (%)",
     col = c("#B2182B"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")],
      type= 'b',
      lwd = 2,
      col = c("#B2182B"),
      pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,
       lwd = 2,
       col =c("#224444"))
dev.off()



#C index
library(pec)
f <- coxph(Surv(`time to event`,event) ~ REEP2 +TMSB15A+DSEL+ID4,data=pbc)###计算出的C-index
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index
plot(c_index)

# 95% confidence interval for c-index
c_index_low <- c_index[1] -c_index[2]
c_index_high <- c_index[1] + c_index[2]




# 计算模型基因和m6A因子的相关系数。

# 进行spearman 相关性分析，返回相关性系数和p值
exprSet_m6A <-  read.table("01rawData/GSE102349_TPM_symbol.txt",header = T,sep = "\t",check.names = F,row.names = 1)
HYPOXIA <- read.table("03Cox/GSE102349_m6A_uniCox.txt",header = T,sep="\t")
m6A_gene <- HYPOXIA$gene[HYPOXIA$coxPvalue < 0.05]
Cluster_res <- read.table("01rawData/GSE102349_Cluster_res.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]
exprSet_m6A <- exprSet_m6A[,c(Cluster1,Cluster2)]
tcga_gsva <- exprSet_m6A[rownames(exprSet_m6A)%in%m6A_gene, ]

tcga_expr <- exprSet_m6A

genelist <- c("REEP2","TMSB15A","DSEL","ID4")
immuscore <- function(gene){
   y <- as.numeric(tcga_expr[gene,])
   colnames <- rownames(tcga_gsva)
   do.call(rbind,lapply(colnames, function(x){
      dd  <- cor.test(as.numeric(tcga_gsva[x,]), y , method="pearson")
      data.frame(gene=gene,m6A_gene=x,cor=dd$estimate,p.value=dd$p.value )
   }))
}


immuscore("REEP2")
# The correlation between genelist and immune infiltration was calculated in batches
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)

write.table(data, "m6A-correlation/m6A10_modelGenes_cor.txt", quote = F, row.names = F,sep = "\t")

data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]

ggplot(data, aes( m6A_gene,gene)) + 
   geom_tile(aes(fill = cor), colour = "white",size=1)+
   scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
   geom_text(aes(label=pstar),col ="black",size = 5)+ #coord_flip()+
   theme_minimal()+
   theme(axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1),
         axis.text.y = element_text(size = 8))+
   
   labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
ggsave("m6A-correlation/m6A_gene10_modelGenes_correlation.pdf", width = 6, height = 2.5)
