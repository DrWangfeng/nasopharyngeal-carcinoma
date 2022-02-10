library(survival)
library(survminer)
inputfile<- "model/validSet_riskScore.txt"
outputfile <- "validSet"
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$event <- ifelse(rt$event =="Last follow-up",0,1)
diff=survdiff(Surv(`time to event`, `event`) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)
HR = (diff$obs[2]/diff$exp[2])/(diff$obs[1]/diff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))

fit <- survfit(Surv(`time to event`, `event`) ~ risk, data = rt)
ggsurvplot(fit,
           pval = TRUE,
           linetype = "solid",  
           palette = mycol[2:3],
           #surv.median.line = "hv", 
           title = "GSE102349",
           ylab = "Progression-free survival (percentage)",
           xlab = " Time (Months)",
           # legend.title = "Survival Plot",
           legend = c(0.8,0.30),
           legend.labs = c("High Risk","Low Risk"),
           legend.title="",
           risk.table = T,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)

dev.copy2pdf(file = paste0("ROC/",outputfile,"_survival.pdf") , width = 5.5,height = 5.5)
dev.off()


library(survivalROC)
period_time <- 12
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$event <- ifelse(rt$event =="Last follow-up",0,1)
rocCol=c("#3cb346","#eeb401", "#ef1828","#942d8d")
aucText <- c()
pdf(file=paste0("ROC/",outputfile,"_ROC.pdf"),width=5,height=5)
par(mar=c(4,4,2,1),mgp=c(2,0.5,0))
roc=survivalROC(Stime=rt$`time to event`,
                status=rt$event,
                marker = rt$riskScore ,
                predict.time =period_time ,
                span = 0.25*nrow(rt)^(-0.20),##span,NNE
                method="NNE")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#1B9E77",
     xlab="1-Specificity (False Positive)", ylab="Sensitivity (True Positive)",
     #main="ROC curve, Method = KM",
     main=outputfile,
     lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1,lty= 3)
j <- 0
for(i in c(1:3)){
        
        roc=survivalROC(Stime=rt$`time to event`, status=rt$event, marker = rt$riskScore,
                        span = 0.25*nrow(rt)^(-0.20),##span,NNE
                        predict.time =period_time*i, method="NNE")
        
        j=j+1
        aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
        lines(roc$FP, roc$TP, type="l",col=rocCol[j],xlim=c(0,1), ylim=c(0,1),lwd = 1.8, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
        
}
legend("bottomright", aucText,lty= 1,lwd=1.8,bty="n",col=rocCol)
dev.off()


## Training set risk curve

library(pheatmap)
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$event <- ifelse(rt$event =="Last follow-up",0,1)
rt=rt[order(rt$riskScore),]                               

riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="Low Risk"])
highLength=length(riskClass[riskClass=="High Risk"])
line=rt[,"riskScore"]

pdf(file= paste0("ROC/",outputfile,"_riskScore.pdf"),width = 6,height = 3)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#028846",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High Risk", "Low Risk"),bty="n",pch=19,col=c("red","#028846"),cex=1.2)
dev.off()


color=as.vector(rt$event)
color[color==1]="red"
color[color==0]="#028846"
pdf(file=paste0("ROC/",outputfile,"_survStat.pdf"),width = 6,height = 3)
plot(rt$`time to event`,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Progression-free time (Months)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","#028846"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()


rt1=rt[c(12:(ncol(rt)-2))]

rt1 <- log2(rt1+1)
rt1=t(rt1)


annotation=data.frame(Type=(rt[,ncol(rt)]))
rownames(annotation)=rownames(rt)
ann_colors = list(Type =c(`Low Risk`="#00BFC4", `High Risk`="#F8766D"))


pdf(file=paste0("ROC/",outputfile,"_heatmap.pdf"),width = 6.5,height = 2.5)
pheatmap(rt1, 
         scale = "row",
         annotation_col=annotation, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         border=FALSE,
         color = colorRampPalette(c("#028846", "white", "red"))(50) )
dev.off()


