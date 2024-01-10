rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/RE/FJM/")


load(file = './Result/COREAD.Rdata')

load(file = './Result/meth001.Rdata')
DEG$Symbols=rownames(DEG)
U_D_gene=unlist(DEG$Symbols)
t1=DEG[,c(7,8)]
t2=pdeg[,c(14,25)]
Hyper_Hypo=unlist(t2$gene)

Hyper_D=intersect(unlist(pdeg[pdeg$change=='Hypermethylation',][14]),unlist(rownames(DEG[DEG$change=='Down-regulated',])))
Hypo_U=intersect(unlist(pdeg[pdeg$change=='Hypomethylation',][14]),unlist(rownames(DEG[DEG$change=='Up-regulated',])))

ACE=c(Hypo_U,Hyper_D)


## load package data.table to speed up reading and saving data.
library(data.table)
library(dplyr)
library(tidyverse)

#read data####

TCGA.COAD.GDC_phenotype <- fread("./data/TCGA-COAD/TCGA-COAD.GDC_phenotype.tsv.gz",header = T, sep = '\t',data.table = F)
TCGA.COAD.htseq_counts <- fread("./data/TCGA-COAD/TCGA-COAD.htseq_counts.tsv.gz",header = T, sep = '\t',data.table = F)
TCGA.COAD.probeMap<- fread("./data/TCGA-COAD/gencode.v22.annotation.gene.probeMap",header = T, sep = '\t',data.table = F)
TCGA.COAD.probeMap=TCGA.COAD.probeMap[,c(1,2)]
head(TCGA.COAD.probeMap)
###merge the transformed probe information with gene expression profile
count=merge(TCGA.COAD.probeMap,TCGA.COAD.htseq_counts,by.x = "id",by.y ="Ensembl_ID" )
rowMan=apply(count[,3:ncol(count)],1,function(x)mean(as.numeric(x),na.rm=T))
count =count[order(rowMan,decreasing =T),]
count =count[!duplicated(count$gene),]
rownames(count)=count$gene
count=count[,3:ncol(count)]
count=round(2**count-1,0)
head(colnames(count))
colnames(count) = gsub("[.]","-",colnames(count))
head(colnames(count))
count=count[,which(substr(colnames(count),14,16) == "01A" | substr(colnames(count),14,16) == "11A" )]#01–09是癌症，10–19是正常，20–29是癌旁

# rownames(count)=as.data.frame(t(as.data.frame(strsplit(rownames(count),"[.]"))))$V1

save(count,ACE,file="./result/modelanddata.Rdata")



a=substr(colnames(count),14,16)
table(a)
remove(a)

#############
coldata=as.data.frame(colnames(count))
colnames(coldata)="sample"
coldata$group=ifelse(substr(coldata$sample,14,16) == "01A","T","N")
mm=model.matrix(~0 + coldata$group)
colnames(mm)=c("N","T")

#DEseq####
library(DESeq2)
dds= DESeqDataSetFromMatrix(countData = count ,colData = coldata,design= ~group)
keep=filterByExpr(count,design=mm )
table(keep)
dds=dds[keep,]
dds=estimateSizeFactors(dds)
norm_count=counts(dds, normalized=TRUE)## normalizing data for down stream analysis
norm_count=log2(norm_count+1)
boxplot(norm_count[,1:30])

###############

library(edgeR)
library(limma)

d=DGEList(counts = count,group = coldata$group)
keep <- filterByExpr(count, design = mm)
table(keep)
d =d[keep,keep.lib.sizes=F] 
d=calcNormFactors(d,method = "TMM")############# optional, may be different with the real data processing
v =voom(counts = d, design =mm, plot=TRUE)
norm=as.data.frame(v$E)##normalizing data for down stream analysis
boxplot(norm[,1:30])


library(edgeR)
dge <- DGEList(counts =count, group = coldata$group)
keep <- rowSums(cpm(dge) > 1 ) >= 2
table(keep)
dge <- dge[keep,keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")############# another optional normalization method
lcpm=cpm(dge,log=T)  ##normalized data for down stream analysis
boxplot(lcpm[,1:30])

###################

surv=read.delim("./data/TCGA-COAD/TCGA-COAD.survival.tsv", row.names=1, stringsAsFactors=FALSE)
# surv$sample=gsub("-",".",rownames(surv))


#### gene expression and survive data

library('glmnet')
library("survminer")
library("survival")
library("openxlsx")
library("GSVA")
set.seed(100)

load( './result/COREAD.Rdata')
remove(COREAD)

norm_count=as.data.frame(norm_count)

norm_count=norm_count[,which(substr(colnames(norm_count),14,16) == "01A")]

COAD=na.omit(norm_count[ACE,])
surv=na.omit(surv[colnames(COAD),])
COAD=COAD[,rownames(surv)]

COAD=as.data.frame(t(COAD))

a=substr(colnames(COAD),14,16)
table(a)
remove(a)

COAD$futime=surv[rownames(COAD),3]
COAD$fustat=surv[rownames(COAD),1]

rt_dexp=COAD

# write.table(rt_dexp,'./Result/rt.txt',sep='\t',row.names = T,col.names = T,quote = F)
# rt=read.table("./Result/rt.txt", header=T, sep="\t", row.names=1, check.names=F)



library(survival)
pFilter=0.05

sigGenes=c("futime","fustat")
outTab=data.frame()

for (i in colnames(rt_dexp[,1:(ncol(rt_dexp)-2)])){            #extract the gene name
  cox<- coxph(Surv(futime,fustat)~rt_dexp[,i], data = rt_dexp)
  coxSummary = summary (cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,cbind(id=i,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
  }
}

uniSigExp=rt_dexp[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
# write.table(uniSigExp, file="./Result/uniSigExp.txt",sep="\t",row.names=F,quote=F)

write.table(outTab, file="./Result/outTabcox.csv",sep=",",row.names=F,quote=F)


if(F){
  
  de_exp=uniSigExp[,c(4,5,6,7,3,2)]
  colnames(de_exp)[5:6]=c('os','os_time')
  de_exp1=de_exp
  p_hr=data.frame()
  numb=ncol(de_exp)-2
  for(i in 1:numb){
    cox_res=summary(coxph(Surv(de_exp$os_time, de_exp$os)~de_exp[,i],data = de_exp))
    p_hr=rbind(p_hr,data.frame(gene=colnames(de_exp)[i],
                               p_rna=as.data.frame(cox_res[10])[3,1],
                               hr_rna=as.data.frame(cox_res[8])[1,1],
                               low_CI=as.data.frame(cox_res[8])[1,3],
                               up_CI=as.data.frame(cox_res[8])[1,4]))
  }
  
  # de_exp1$OS_months=as.numeric(de_exp1$OS_months)
  # de_exp1$status_n=as.numeric(de_exp1$status_n)
  
  
  for(i in 1:numb){de_exp1[,i]=ifelse(de_exp1[,i] > median(de_exp1[,i]),"high","low")}
  p_km=list()
  for(i in 1:numb){p_km=append(p_km,surv_pvalue(survfit(Surv(de_exp1$os_time, de_exp1$os)~de_exp1[,i],data = de_exp1))$pval)}
  p_km=data.frame(gene=colnames(de_exp)[1:numb],km=unlist(p_km))
  
  de_km_cox=merge(p_km,p_hr,by="gene")
  
  colnames(de_km_cox)=c("gene","km_p","cox_p","cox_hr","hr_low","hr_high")
  rownames(de_km_cox)=de_km_cox$gene
  library("survival")
  library("survminer")
  zz=de_km_cox[de_km_cox$km_p<0.05,'gene']
  ggsurvplot(survfit(as.formula(paste('Surv(os_time, os) ~',zz[3])), data = de_exp1),pval = TRUE, conf.int = FALSE,risk.table = TRUE,risk.table.col = "strata",linetype = "strata",surv.median.line = "hv",ggtheme = theme_bw(),palette = c("#E7B800", "#2E9FDF"))
  
}





library("glmnet")
library ("survival")
library(survminer)

rt_uniSigExp=uniSigExp[,-1]
rt_uniSigExp$futime[rt_uniSigExp$futime<=0]=1  # when the follow up time is 0, code will report error, so all the zeros changed into 1
x=as.matrix(rt_uniSigExp[,c(3:ncol(rt_uniSigExp))])
y=data.matrix(Surv(rt_uniSigExp$futime,rt_uniSigExp$fustat))


fit<-glmnet(x,y,family = "cox", maxit = 1000)#lasso regression  randomly simulation of 1000 times.
plot(fit,xvar = "lambda", label = TRUE)

cvfit<- cv.glmnet(x,y,family="cox",maxit = 1000)# cross-validation
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.lse)),lty='dashed')# mark with the little dot


coef <-coef(fit,s = cvfit$lambda.min)
index <- which(coef[,1]!= 0)# extract the coefficients with non-zeros ##### source codes as index <- which(coef!= 0)
actCoef <- coef[index]
lassoGene=row.names (coef)[index]
lassoGene=c("futime","fustat", lassoGene)
lassoSigExp=rt_uniSigExp[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
# write.table(lassoSigExp,file="./Result/lassoSigExp.txt",sep="\t", row.names=F, quote=F)

rt_lassoSigExp=lassoSigExp[,-1]

rt_lassoSigExp[,"futime"]=rt_lassoSigExp[,"futime"] /31
#使用train组构建COX模型  多因素
multiCox=coxph(Surv(futime,fustat)~.,data = rt_lassoSigExp)
# multiCox=step(multiCox, direction = "both")# select gene both forward and backward, lasso depend on correlation, both direction results in 5 genes
# multiCox=step(multiCox, direction = "backward")
multiCox=step(multiCox, direction = "forward")# forward gene selection, lasso accords to correlation, forward results in 7 genes
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])


outTab=cbind(id=row.names(outTab),outTab)

ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22,0.4),
         fontsize = 0.7,
         refLabel = "reference",noDigits = 2)


riskScore=predict(multiCox,type="risk",newdata=rt_lassoSigExp)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub ("'","",coxGene)# pay attention to the "、" marker
outCol=c("futime","fustat",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))


rt=cbind(id=rownames(cbind(rt_lassoSigExp[,outCol],riskScore,risk)),cbind(rt_lassoSigExp[, outCol],riskScore,risk))

rt_riskTrain=rt[,-1]



diff=survdiff(Surv(futime,fustat)~risk,data = rt)
pValue=1-pchisq(diff$chisq, df=1)
pValue=signif (pValue,4)
pValue=format (pValue,scientific = TRUE)
fit <- survfit (Surv(futime,fustat)~risk, data = rt)
summary(fit) # analyze the five-year survival 
plot (fit,
      lwd=2,
      col=c ("red","blue"),
      xlab="Time (year)",ylab="Survival rate",
      main=paste("Survival curve (p=",pValue,")",sep=""),mark.time=T)

legend('topright',c('high risk','low risk'),lwd=2,col=c('red','blue'))


library(survivalROC)
par(oma=c(0.5,1,0,1) , font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt_riskTrain$futime, status=rt_riskTrain$fustat,marker = rt_riskTrain$riskScore,
                predict.time =3,method="KM")
plot(roc$FP,roc$TP, type="l",xlim=c(0,1),ylim=c(0,1),col='red',xlab="False positive rate", ylab="True positive rate",
     main=paste("RoC curve (","AUC = ",round(roc$AUC,3),")"),lwd = 2,cex.main=1.3,cex.lab=1.2,cex.axis=1.2,font=1.2)

abline(0,1)

