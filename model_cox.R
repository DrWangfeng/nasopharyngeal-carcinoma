rm(list=ls())

#install.packages("survival")

library(survival)
library(tidyverse)
clincal <- read.table("01rawData/GSE102349_clinical_PFS.txt",sep = "\t",header = T,check.names = F)

expFile="01rawData/GSE102349_TPM_symbol.txt" 
expr=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)

DEGs_File="DEGs/GSE102349_edgeR_diff.txt" 
DEGs_File=read.table(DEGs_File,sep="\t",header=T,check.names=F,row.names = 1)

exprSet_m6A <- expr[rownames(DEGs_File),]

exprSet_m6A1 <- as.data.frame( t(exprSet_m6A))

exprSet_m6A1$id <- substr(rownames(exprSet_m6A1),1,12)
exprSet_m6A1_Time <- merge(clincal,exprSet_m6A1,by="id")

exprSet_m6A1_Time[1:3,1:4]
write.table(exprSet_m6A1_Time,file = "model/GSE102349_DEGs_exprSetTime.txt",row.names = F,quote = F,sep = "\t")

# Univariate Cox analysis
rt=read.table("model/GSE102349_DEGs_exprSetTime.txt",header=T,check.names=F,row.names=1,sep = "\t") 
whole <- rt1 <- rt

rt$event <- ifelse(rt$event =="Last follow-up",0,1)

# Select the name of the variable to be used for univariate analysis
(variable_names<-colnames(rt)[-c(1:10)])


library("caret")
# Set the seed so that it can be repeated
set.seed(359)
# 7:3 Randomly divided into training set and test set
trainId  <- createDataPartition(y=rt$`time to event`,p=0.7,list = FALSE)

trainSet <- rt[trainId,]
testSet <- rt1[-trainId,]


rt <- trainSet
rt1 <- rt1[trainId,]

library(survival)
sur<-Surv(time=rt$`time to event`, event = rt$event)

pFilter=0.05
outTab=data.frame()
sigGenes = colnames(rt1)[1:10]
for(gene in variable_names){
  if(sd(rt[,gene])<1){next}
  if(grepl("-", gene)){next}
  cox=coxph(as.formula(paste0('sur~',gene))  ,data = rt)
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
                         KM=pValue,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         coxPvalue=coxP) )
    }
  }
}

write.table(outTab,file="model/trianSet_uniCox.txt",sep="\t",row.names=F,quote=F)   

surSigExp=rt1[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="model/trianSet_uniSigExp.txt",sep="\t",row.names=F,quote=F)

# Draw the forest map function
bioForest=function(coxFile=null,forestFile=null){
  # Reading input file
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <-  gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\", rownames(rt)) 
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))
  
  # Output graphics
  pdf(file=forestFile, width = 5,height = 8)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(5,2))
  
  # Plot the clinical information on the left side of the forest map
  xlim = c(0,3)
  par(mar=c(4,3,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.6-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.6-0.5*0.2,n+1,'Pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard Ratio',cex=text.cex,font=2,adj=1,)
  

  par(mar=c(4,0,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)+0.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard Ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=1.5)
  abline(v=1,col="black",lty=2,lwd=1.3)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.8)
  axis(1)
  dev.off()
}

bioForest(coxFile="model/trianSet_uniCox.txt",forestFile="model/trian_uniForest.pdf")

####################################
# Multivariate Cox analysis
options(stringAsFactors=F)
library(survival)  
library(glmnet)

rt=read.table("model/trianSet_uniSigExp.txt",header=T,sep="\t",check.names=F,row.names = 1)    #读取输入文件
rt1 = rt
rt$event <- ifelse(rt$event =="Last follow-up",0,1)

rt <- rt[rt$`time to event` != 0,]

str(rt)
# Multivariate Cox +LASSO penalty
set.seed(9)

x = as.matrix(rt[,-c(1:10)])
y = Surv(time=rt$`time to event`, event = rt$event)
cvfit = cv.glmnet(x,y, family = "cox", nfold = 10)
pdf("model/cvfit.pdf",width = 6,height = 6)
plot(cvfit)
# abline(v = c(log(cvfit$lambda.min), log(cvfit$lambda.1se)),lty=2)+
text(x = log(cvfit$lambda.min),y = 20,
     paste('Lambda.min\n',round(cvfit$lambda.min,4)),cex=1,adj=0.9)
text(x = log(cvfit$lambda.1se)+0.01,y = 40,
     paste('Lambda.lse\n',round(cvfit$lambda.1se,4)),cex=1,adj=0.9)
dev.off()

fit <- glmnet(x, y, family = "cox",nfold = 10)
pdf(file="model/fit.pdf",width=6,height=5.5)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()


myCoefs <- coef(cvfit, s=cvfit$lambda.min)
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea
sigGenes= colnames(rt1)[1:10]
rt <- rt1[,c(sigGenes,lasso_fea)]
rt$riskScore <-  apply(rt[,lasso_fea], 1, function(x) {x %*% myCoefs@x})


rt$risk=as.vector(ifelse(rt$riskScore>median(rt$riskScore),"High Risk","Low Risk"))

write.table(cbind(id=rownames(rt),rt),
            file="model/trianSet_riskScore.txt",sep="\t",quote=F,row.names=F)

lassoGene <-  cbind(Gene = lasso_fea,Coef = myCoefs[which(myCoefs != 0 )])
write.table(lassoGene,file="model/trianSet_lasso_coef.txt",row.names = F,quote = F,sep = "\t")

testSet <- testSet[,c(sigGenes,lasso_fea)]
testSet$riskScore <- apply(testSet[,lasso_fea], 1, function(x) {x %*% myCoefs@x})
testSet$risk=as.vector(ifelse(testSet$riskScore>median(testSet$riskScore),"High Risk","Low Risk"))
write.table( cbind(id=rownames(testSet),testSet),
             file="model/testSet_riskScore.txt",
             sep="\t",row.names = F,
             quote=F)

# total data set
whole <- whole[,c(sigGenes,lasso_fea)]
whole$riskScore <- apply(whole[,lasso_fea], 1, function(x) {x %*% myCoefs@x})

whole$risk=as.vector(ifelse(whole$riskScore>median(whole$riskScore),"High Risk","Low Risk"))

write.table( cbind(id=rownames(whole),whole),
             file="model/wholeSet_riskScore.txt",
             sep="\t",row.names = F,
             quote=F)
