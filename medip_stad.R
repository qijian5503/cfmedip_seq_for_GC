rm(list=ls())
gc()
options(stringsAsFactors = F)
setwd("~/MEDIP/medip_stad/")

library(limma)
library(edgeR)
library(ROCR)
library(pROC)
library(caret)
library(glmnet)
library(dplyr)
library(randomForest)
library(Rtsne)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(caret)
library(gtsummary)
library(tableone)
library(MatchIt)

##########################################################################
## peak data           ##########################
peak_exp <- read.delim("peaks_exp.bed", header=FALSE)
peak_exp$V9=toupper(gsub("-","_",peak_exp$V9))
peak_exp$peak=paste(peak_exp$V1,peak_exp$V2,peak_exp$V3,sep = "_")
peak_exp$a=paste(peak_exp$peak,peak$V9,sep = "_")
peak_exp=peak_exp[order(peak_exp$V8,decreasing = T),]
peak_exp=peak_exp[!duplicated(peak_exp$a),]

peak_exp=reshape2::dcast(peak_exp,peak~V9,value.var = "V8")
rownames(peak_exp)=peak_exp$peak
peak_exp$peak=NULL
peak_exp[is.na(peak_exp)]=0

sel_peak=str_split(rownames(peak_exp),"_",simplify = T) %>% as.data.frame()
sel_peak=sel_peak[sel_peak$V4 %in% "" & !sel_peak$V1 %in% c("chrM","chrUn","chrX","chrY"),]
peak_exp=peak_exp[rownames(peak_exp) %in% paste(sel_peak$V1,sel_peak$V2,sel_peak$V3,sep = "_"),]
#save(peak_exp,file="medip_stad_peak_samp250_raw_exp.RData")


load("medip_stad_peak_samp250_phen.RData")
train_test=phen[phen$group %in% "train_test",]
train_test_samp=train_test$sample
subsample=createDataPartition(train_test$Group, p = 0.8,times = 100)
for (i in 1:100) {
    train_samp=train_test_samp[subsample[[i]]]
    test_samp=train_test_samp[!train_test_samp %in% train_samp]
    subsample[[i]]=list(train_samp=train_samp,test_samp=test_samp)
}
#save(subsample,file="medip_stad_peak_samp250_subsamp.RData")

## peak annotation     ##########################################
rm(list=ls())
options(stringsAsFactors = F)

load("medip_stad_peak_samp250_raw_exp.RData")
peak_ana=data.frame(pos=rownames(peak_exp),max=rowSums(peak_exp),count=rowSums(peak_exp>0))
peak_ana=peak_ana[!grepl("chrUn|random|hap",peak_ana$pos),]

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(tidyverse)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene

pos=do.call(rbind,lapply(as.character(peak_ana$pos), function(x){unlist(strsplit(x,"_"))})) %>% as.data.frame()
pos=na.omit(data.frame(seqnames=pos$V1,start=as.numeric(pos$V2),end=as.numeric(pos$V3)))
peakAnno = GRanges(seqnames=Rle(pos[,1]),ranges=IRanges(pos[,2], pos[,3]), strand=rep(c("*"), nrow(pos))) %>% annotatePeak(. ,TxDb=txdb, annoDb="org.Hs.eg.db") %>% as.data.frame()
peakAnno$annotation1=str_split(peakAnno$annotation,'[(]',simplify = T)[,1]
table(peakAnno$annotation1)

peakAnno$pos=paste(peakAnno$seqnames,peakAnno$start,peakAnno$end,sep = "_")
peak_ana=merge(peak_ana,peakAnno[,c(6,14:16,18:19)],by="pos")
#save(peak_ana,file="medip_stad_peak_samp250_peakAnno.RData")

## sample QC           ##########################################
rm(list=ls())
options(stringsAsFactors = F)

load("medip_stad_peak_samp250_all_coverage.RData")
load("medip_stad_peak_samp250_raw_exp.RData")
load("medip_stad_peak_valid_subsamp.RData")
load("medip_stad_peak_samp250_phen.RData")

train_test=phen[phen$group %in% "train_test",]
train_test_samp=train_test$sample
table(train_test_samp %in% colnames(peak_exp))
table(phen$sample %in% colnames(peak_exp))

peak_tab=as.data.frame(colSums(peak_exp>0))
summary(peak_tab$`colSums(peak_exp > 0)`)
rownames(all_coverage)=all_coverage$V1
table(rownames(peak_tab) %in% rownames(all_coverage))
peak_tab=cbind(peak_tab[all_coverage$V1,],all_coverage)
colnames(peak_tab)[1]="peak_tab"
peak_tab=cbind(phen,peak_tab[rownames(phen),])
samp_select=rownames(peak_tab[peak_tab$peak_tab> 1000 & peak_tab$V2>10000000,])
#save(samp_select,file = "medip_stad_peak_samp250_samp15_select.RData")


## peak norm           ##########################################
rm(list=ls())
options(stringsAsFactors = F)
library(limma)
library(edgeR)

load("medip_stad_peak_samp250_raw_exp.RData")
load("medip_stad_peak_samp250_samp15_select.RData")
load("medip_stad_peak_samp250_phen.RData")
table(phen$group)
train_test=phen[phen$group %in% "train_test" & phen$sample %in% samp_select,]
train_test_samp=train_test$sample

train_test=peak_exp[,train_test_samp]
keep <- rowSums(train_test > 0) >= ncol(train_test)*0.2
table(keep)
#save(keep,file = "medip_stad_peak_samp250_keep_10_02.RData")

peak=rownames(train_test)[keep]
peak=str_split(peak,"_",simplify = T)
#write.table(peak,"peak_keep_name.bed",quote = F,row.names = F,col.names = F,sep = "\t")

train_test=train_test[keep,]
f75 <- rep_len(1,ncol(train_test))
for (j in seq_len(ncol(train_test))) f75[j] <- quantile(train_test[,j], probs=0.75)
lib.size=colSums(train_test)
f75=f75 / lib.size
refColumn <- colnames(train_test)[which.min(abs(f75-mean(f75)))]
#save(refColumn,file="medip_stad_peak_samp250_refColumn_10_02.RData")

ref_exp=data.frame(peak=rownames(peak_exp),refColumn=peak_exp[,refColumn])
#save(ref_exp,file = "ref_exp.RData")

d=DGEList(counts = train_test) %>% calcNormFactors(. ,method = "TMM" ,refColumn = refColumn) %>% voom()
norm=as.data.frame(d$E)
norm[1:3,1:3]
#save(norm,file = "medip_stad_peak_samp250_traintest_norm_voom_10_02_noxy.RData")

## validation norm     #####################################################
valid_peak_exp <- read.delim("valid_peaks_exp.bed", header=FALSE)
valid_peak_exp$V9=toupper(gsub("-","_",valid_peak_exp$V9))
valid_peak_exp$peak=paste(valid_peak_exp$V1,valid_peak_exp$V2,valid_peak_exp$V3,sep = "_")
valid_peak_exp$a=paste(valid_peak_exp$peak,peak$V9,sep = "_")
valid_peak_exp=valid_peak_exp[order(valid_peak_exp$V8,decreasing = T),]
valid_peak_exp=valid_peak_exp[!duplicated(valid_peak_exp$a),]
valid_peak_exp=reshape2::dcast(valid_peak_exp,peak~V9,value.var = "V8")
valid_peak_exp[is.na(valid_peak_exp)]=0

load("ref_exp.RData")
valid=merge(valid_peak_exp,ref_exp,by="peak",all.y=T)
rownames(valid)=valid$peak
valid$peak=NULL
valid[is.na(valid)]=0

d=DGEList(counts = valid) %>% calcNormFactors(. ,method = "TMM" ,refColumn = "V2") %>% voom()
norm_v=as.data.frame(d$E)
table(norm_v[,"refColumn"]==norm[,refColumn])
norm_v$refColumn=NULL
#save(norm_v,file = "medip_stad_peak_samp250_norm_voom_10_02_valid_noxy.RData")


## DMR                 #################################################################
rm(list=ls())
options(stringsAsFactors = F)

load("medip_stad_peak_samp250_subsamp.RData")
load("medip_stad_peak_samp250_traintest_norm_voom_10_02_noxy.RData")
load("medip_stad_peak_samp250_samp15_select.RData")

load("medip_stad_peak_samp250_phen.RData")
phen=na.omit(phen[samp_select,c(2,3,5)])
phen$Gender=as.numeric(phen$Gender)-1

limma_de_all=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(phen)]
    coldata=phen[train_samp,]
    psm <- matchit(Group~Age+Gender,data=coldata,method="nearest",distance = "logit",ratio = 1)
    coldata <- match.data(psm)
    design=model.matrix(~0+coldata$subclass+coldata$Group)
    colnames(design)=gsub("coldata[$]","",colnames(design))
    res=lmFit(norm[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="Groupcancer") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[abs(res$logFC)>3 & res$adj.P.Val<0.0001,]
    message(nrow(limma_de_all[[i]]))
}
#save(limma_de_all,file="medip_stad_peak_samp250_norm_voom_10_02_limma_PSM_de_noxy_all.RData")

## lasso               ######################################################
rm(list=ls())
load("medip_stad_peak_samp250_peakAnno.RData")
load("medip_stad_peak_samp250_norm_voom_10_02_limma_PSM_de_noxy_all.RData")

de_gene=data.frame()
for(i in 1:100){
    res=limma_de_all[[i]]
    res$pos=rownames(res)
    de_gene=rbind(de_gene,res)
}
rownames(de_gene)=NULL
de_gene$group=ifelse(de_gene$logFC>0,"Hypermethylation","Hypomethylation")
de_gene=merge(de_gene,peak_ana[,c(1,7,8)],by="pos",all.x=T)

de_pos=as.data.frame(table(de_gene$pos,de_gene$group))
de_pos$Var1=as.character(de_pos$Var1)
de_pos=reshape2::dcast(de_pos,Var1~Var2,value.var = "Freq")
table(de_pos$Hypermethylation>0 , de_pos$Hypomethylation>0)
de_pos=de_pos[de_pos$Hypermethylation>10 | de_pos$Hypomethylation>10,]
de_pos=as.character(de_pos$Var1)

library(glmnet)
de_exp=as.data.frame(t(norm[de_pos,]))
de_exp$Group=phen[rownames(de_exp),"Group"]
fit <- glmnet(as.matrix(de_exp[,de_pos]), de_exp$Group, family = "binomial",type.measure = "class")
#save(fit,file = "lasso_fit.RData")

cv_fit = cv.glmnet(as.matrix(de_exp[,de_pos]), de_exp$Group, family = "binomial",type.measure = "class")
coef.min = coef(cv_fit, s = "lambda.min")
gen=coef.min@Dimnames[[1]][coef.min@i+1][-1] 
#save(cv_fit,file = "lasso_gene21.RData")



## model               #####################################################
rm(list=ls())

load("medip_stad_peak_samp250_subsamp.RData")
load("medip_stad_peak_samp250_phen.RData")
load("medip_stad_peak_samp250_samp15_select.RData")
load("medip_stad_peak_samp250_traintest_norm_voom_10_02_noxy.RData")
load("medip_stad_peak_samp250_norm_voom_10_02_valid_noxy.RData")
valid_sample=colnames(norm_v)

load("lasso_gene21.RData")
coef.min = coef(cv_fit, s = "lambda.min")
gen=coef.min@Dimnames[[1]][coef.min@i+1][-1] 

library(glmnet)
library(randomForest)
ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% samp_select]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% samp_select]
    de_exp=as.data.frame(t(norm[gen,]))
    de_exp$group=factor(phen[rownames(de_exp),"Group"],levels = c("cancer","normal"))
    train_data=de_exp[train_samp,c(gen,"group")]
    test_data=de_exp[test_samp,c(gen,"group")]
    
    rf_mod = train(group ~ .,data = train_data,method = "rf",trControl = ctrl)
    rf_pre=predict(rf_mod,test_data,type="prob")%>%data.frame
    rf_auc=unlist(slot(performance(prediction(rf_pre$normal,test_data$group),'auc'),"y.values"))
    valid_data=as.data.frame(t(norm_v[gen,]))
    valid_pre=predict(rf_mod,valid_data,type="prob")%>%data.frame
    valid_pre$group=as.factor(phen[rownames(valid_pre),"Group"])
    valid_auc=unlist(slot(performance(prediction(valid_pre$normal,valid_pre$group),'auc'),"y.values"))
    
    limma_pre=rbind(data.frame(sample=rownames(test_data),group=test_data$group,pre=rf_pre$normal),data.frame(sample=rownames(valid_pre),group=valid_pre$group,pre=valid_pre$normal))
    limma_modlist[[i]]= list(limma_pre=limma_pre,rf_mod=rf_mod,limma_de=de,valid_auc=valid_auc,rf_auc=rf_auc)
    message(round(rf_auc,3),"   ",round(valid_auc,3))
}
#save(limma_modlist,file="medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21.RData")


## model with cfdna    #####################################################
rm(list=ls())

load("medip_stad_peak_samp250_subsamp.RData")
load("medip_stad_peak_samp250_phen.RData")
load("medip_stad_peak_samp250_samp15_select.RData")
load("medip_stad_peak_samp250_traintest_norm_voom_10_02_noxy.RData")
load("medip_stad_peak_samp250_norm_voom_10_02_valid_noxy.RData")
valid_sample=colnames(norm_v)


load("lasso_gene21.RData")
coef.min = coef(cv_fit, s = "lambda.min")
de_gene=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因

de_exp=as.data.frame(t(norm[de_gene,]))
de_exp$group=factor(phen[rownames(de_exp),"Group"],levels = c("cancer","normal"))
de_exp$Concentration=phen[rownames(de_exp),"Concentration"]
de_exp=na.omit(de_exp)

valid_data=as.data.frame(t(norm_v[de_gene,]))
valid_data$Concentration=phen[rownames(valid_data),"Concentration"]
valid_data=na.omit(valid_data)

library(glmnet)
library(randomForest)
ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(de_exp)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(de_exp)]
    train_data=de_exp[train_samp,c("Concentration",de_gene,"group")]
    test_data=de_exp[test_samp,c("Concentration",de_gene,"group")]
    
    rf_mod = train(group ~ .,data = train_data,method = "rf",trControl = ctrl)
    rf_pre=predict(rf_mod,test_data,type="prob")%>%data.frame
    rf_auc=unlist(slot(performance(prediction(rf_pre$normal,test_data$group),'auc'),"y.values"))
    
    valid_pre=predict(rf_mod,valid_data,type="prob")%>%data.frame
    valid_pre$group=as.factor(phen[rownames(valid_pre),"Group"])
    valid_auc=unlist(slot(performance(prediction(valid_pre$normal,valid_pre$group),'auc'),"y.values"))
    
    limma_pre=rbind(data.frame(sample=rownames(test_data),group=test_data$group,pre=rf_pre$normal),data.frame(sample=rownames(valid_pre),group=valid_pre$group,pre=valid_pre$normal))
    limma_modlist[[i]]= list(limma_pre=limma_pre,rf_mod=rf_mod,limma_de=de,valid_auc=valid_auc,rf_auc=rf_auc)
    message(round(rf_auc,3),"   ",round(valid_auc,3))
}
#save(limma_modlist,file="medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21_with_cfdna.RData")

#################################################################################################
## plot                #################################################################################################
## phen statistic      ########################################
rm(list=ls())
options(stringsAsFactors = F)

load("medip_stad_peak_samp250_phen.RData")
colnames(phen)=c("Sample","Age","Gender","cfDNA Concentration(ng/ml)","Group","Stage","group2")
unique(phen$Stage)

library(ggpubr)
phen$group2=factor(phen$group2,labels = c("discovery set","validation set"),levels = c("train_test","valid"))
phen$Group=factor(phen$Group,labels = c("Normal","Cancer"),levels = c("normal","cancer"))
phen$Stage=factor(phen$Stage,labels = c("Normal","Stage 0","Stage I","Stage II","Stage III","Stage IV"),levels = c("normal","0","I","II","III","IV"))

load("medip_stad_peak_samp250_all_coverage.RData")
colnames(all_coverage)=c("Sample","Reads Number(Millions)")
rownames(all_coverage)=all_coverage$Sample

load("medip_stad_peak_samp250_raw_exp.RData")
peak_tab=as.data.frame(colSums(peak_exp>0))
summary(peak_tab$`colSums(peak_exp > 0)`)
rownames(all_coverage)=all_coverage$V1
table(rownames(peak_tab) %in% rownames(all_coverage))
colnames(peak_tab)[1]="Peak Number(Thousands)"
peak_tab$Sample=rownames(peak_tab)
peak_tab=merge(peak_tab,all_coverage,by="Sample")
peak_tab=merge(phen,peak_tab,by="Sample")
colnames(peak_tab)
colnames(peak_tab)[2]="Age(years)"
peak_tab$`Peak Number(Thousands)`=round(peak_tab$`Peak Number(Thousands)`/1000,1)
peak_tab$`Reads Number(Millions)`=round(peak_tab$`Reads Number(Millions)`/1000000,1)

load("mod_all_pre.RData")
all_pre_medina=reshape2::dcast(all_pre,sample+group~group,fun.aggregate = median)
all_pre_medina$Risk_Score=ifelse(all_pre_medina$group %in% "normal",all_pre_medina$normal,all_pre_medina$cancer)
all_pre_medina=all_pre_medina[order(all_pre_medina$Risk_Score,decreasing = T),]
rownames(all_pre_medina)=all_pre_medina$sample

peak_tab$Risk_Score=all_pre_medina[match(peak_tab$Sample,all_pre_medina$sample),"Risk_Score"]
#save(peak_tab,file = "peak_tab.RData")

colnames(peak_tab)
fig_1a=ggplot(peak_tab[!is.na(peak_tab$Group),], aes(Group, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=Group))+xlab(label = "")+
    scale_fill_manual(values=c("#0571B0","#CA0020"))+
    stat_compare_means(aes(group=Group),hide.ns = T)+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

sfig_1a=ggplot(peak_tab[!is.na(peak_tab$Gender) & !is.na(peak_tab$Group),], aes(Group, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=Gender))+xlab(label = "")+
    stat_compare_means(aes(group=Gender),hide.ns = T)+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
    scale_fill_manual(values=c("#0571B0","#CA0020","#ffad60"))

sfig_1b=ggplot(peak_tab[!is.na(peak_tab$`Age(years)`) & !is.na(peak_tab$Group),], aes(Group, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=`Age(years)`))+xlab(label = "")+
    stat_compare_means(aes(group=`Age(years)`),hide.ns = T)+labs(title = "")+ylab("Age(years)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
    scale_fill_manual(values=c("#0571B0","#CA0020","#ffad60"))

my_comparisons=list(c("Stage 0","Stage I"),c("Stage 0","Stage II"),c("Stage 0","Stage III"),c("Stage II","Stage III"),c("Stage III","Stage IV"),c("Normal","Stage III"))
fig_1b=ggplot(peak_tab[!is.na(peak_tab$Stage) ,], aes(Stage, `cfDNA Concentration(ng/ml)`)) + 
    geom_boxplot(aes(fill=Stage))+xlab(label = "")+
    stat_compare_means(comparisons = my_comparisons,aes(group=Stage),hide.ns = T)+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))


colnames(peak_tab)
library(gtsummary)
t1=tbl_summary(peak_tab[peak_tab$group2 %in% "discovery set",c(2:6,8,9)],by = Group,
               statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
               digits = all_continuous() ~ 2,#确定小数点数
               missing_text = "(Missing)")
t2=tbl_summary(peak_tab[peak_tab$group2 %in% "validation set",c(2:6,8,9)],by = Group,
               statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
               digits = all_continuous() ~ 2,#确定小数点数
               missing_text = "(Missing)")
tbl_merge(tbls = list(t1, t2),tab_spanner = c("****", "**validation set**"))


## pos annotation      ###############################
rm(list=ls())
load("medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_all_rf_noxy.RData")
load("medip_stad_peak_samp250_peakAnno.RData")
load("medip_stad_peak_samp250_norm_voom_10_02_limma_PSM_de_noxy_all.RData")

imp_de_gene=data.frame()
for(i in 1:100){
    imp=limma_modlist[[i]]$rf_mod$finalModel$importance %>% as.data.frame()
    imp$pos=rownames(imp)
    de=limma_de_all[[i]]
    imp_de_gene=rbind(imp_de_gene,cbind(imp,de[,c(1,5)]))
}
rownames(imp_de_gene)=NULL

imp_de_gene$group=ifelse(imp_de_gene$logFC>0,"Hypermethylation","Hypomethylation")
imp_de_gene=merge(imp_de_gene,peak_ana[,c(1,7,8)],by="pos",all.x=T)
#save(imp_de_gene,file = "imp_de_gene.RData")

load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
imp_de_gene2=imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,]
imp_de_gene2$group=imp_de_gene[match(imp_de_gene2$pos,imp_de_gene$pos),"group"]

hgnc_complete_set <- read.delim("~/database/hgnc_complete_set.txt")
table(imp_de_gene2$SYMBOL %in% hgnc_complete_set$symbol)
imp_de_gene2$SYMBOL[!imp_de_gene2$SYMBOL %in% hgnc_complete_set$symbol ]
imp_de_gene2$type=hgnc_complete_set[match(imp_de_gene2$SYMBOL,hgnc_complete_set$symbol),"locus_group"]
# write.table(imp_de_gene2,"imp_de_gene2.tsv",quote = F,sep = "\t",row.names = F)
# save(imp_de_gene2,file = "imp_de_gene2.RData")


a= imp_de_gene2[,c("type","group")] %>% dplyr::group_by(group,type) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
a$n_lab=ifelse(a$Freq<1,"",a$n)
a$Freq_lab=ifelse(a$Freq<1,"",paste(round(a$Freq,2),"%"))

fig_1c=ggplot(a,aes(x=group,y=Freq,fill=type))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))

b= imp_de_gene2[,c("group","annotation1")] %>% dplyr::group_by(group,annotation1) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
b$n_lab=ifelse(b$Freq<1,"",b$n)
b$Freq_lab=ifelse(b$Freq<1,"",paste(round(b$Freq,2),"%"))

fig_1d=ggplot(b,aes(x=group,y=Freq,fill=annotation1))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))

## Rtsne               ###########################################################
load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
table(imp_de_gene2$Hypermethylation>10)
table(imp_de_gene2$Hypomethylation>10)
de_pos=as.character(imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,"pos"])

load("medip_stad_peak_samp250_traintest_norm_voom_10_02_noxy.RData")
load("medip_stad_peak_samp250_norm_voom_10_02_valid_noxy.RData")
norm=cbind(norm,norm_v)

load("medip_stad_peak_samp250_phen.RData")
library(Rtsne)
dat=as.data.frame(t(norm[de_pos,]))
tsne_out=Rtsne(as.matrix(dat),pca=FALSE)
library(ggplot2)
dat_plt=as.data.frame(tsne_out$Y)
rownames(dat_plt)=rownames(dat)
dat_plt$txt=rownames(dat)
dat_plt$Group=phen[rownames(dat_plt),"Group"]
dat_plt$Group=factor(dat_plt$Group,labels = c("Normal","Cancer"))

fig_1f=ggplot(dat_plt,aes(x=dat_plt[,1],y=dat_plt[,2],color=Group))+geom_point()+
    scale_color_manual(values=c("#0571B0","#CA0020"))+
    theme(legend.title = element_text(size = 13),legend.text = element_text(size = 13),axis.title.y= element_text(size=15),axis.title.x= element_text(size=15),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+
    labs(x = "",y = "",title ="" )



## pheatmap            ###################################################################################
library(pheatmap)
library(gplots)
library(colorRamps)
library(RColorBrewer)

load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
table(imp_de_gene2$Hypermethylation>10)
table(imp_de_gene2$Hypomethylation>10)
de_pos=as.character(imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,"pos"])

load("medip_stad_peak_samp250_traintest_norm_voom_10_02_noxy.RData")
load("medip_stad_peak_samp250_norm_voom_10_02_valid_noxy.RData")
norm=cbind(norm,norm_v)
nmf_exp=norm[de_pos,]

load("medip_stad_peak_samp250_phen.RData")
table(phen$sample %in% colnames(nmf_exp))
phen=phen[phen$sample %in% colnames(nmf_exp),]
colnames(phen)
phen=phen[,c(5,6,3,2)]
phen=phen[order(phen$Group),]

phen$Group=factor(phen$Group,labels = c("Normal","Cancer"),levels = c("normal","cancer"))
phen$Stage=factor(phen$Stage,labels = c("Normal","Stage 0","Stage I","Stage II","Stage III","Stage IV"),levels = c("normal","0","I","II","III","IV"))
median(phen$Age,na.rm = T)
phen$Age=ifelse(phen$Age>=56,">=56",ifelse(phen$Age<56,"<56",NA))
phen$Age=factor(phen$Age,levels = c("<56",">=56"))


ann_colors = list(
    Group=c(Cancer="#CA0020", Normal="#0571B0"),
    change=c(Hypermethylation="#CA0020", Hypomethylation="#0571B0"),
    Gender=c(`F`="#FB9A99",`M`="#A65628"),
    Age=c(`>=56`="#4DAF4A",`<56`="#984EA3"),
    Stage=c(`Normal`="#0571B0",`Stage 0`="#F5F5F5",`Stage I`="#C7EAE5",`Stage II`="#80CDC1",`Stage III`="#35978F",`Stage IV`="#01665E",`Stage V`="#003C30")
)
phen1=data.frame(row.names = rownames(phen),Group=phen$Group)


split_row=unique(imp_de_gene[c(1,5)])
rownames(split_row)=split_row$pos
split_row=split_row[de_pos,]
colnames(split_row)[2]="change"
split_row$pos=NULL

fig_1e=ComplexHeatmap::pheatmap(as.matrix(nmf_exp[,rownames(phen)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors,
                            col = colorRampPalette(c("navy","white","firebrick3"))(10),
                            annotation_col = phen,annotation_row = split_row,
                            row_split = split_row$change,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)





## pos top1000         ######################################
rm(list=ls())
load("medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_all_rf_noxy.RData")
load("medip_stad_peak_samp250_peakAnno.RData")
load("medip_stad_peak_samp250_norm_voom_10_02_limma_PSM_de_noxy_all.RData")
limma_de=list()
for(i in 1:100){
    res=limma_de_all[[i]][1:1000,]
    limma_de=rbind(limma_de,data.frame(res,pos=rownames(res)))
}
rownames(limma_de)=NULL

limma_de$group=ifelse(limma_de$logFC>0,"Hypermethylation","Hypomethylation")
limma_de=merge(limma_de,peak_ana[,c(1,7,8)],by="pos",all.x=T)
#save(limma_de,file = "limma_de_top1000.RData")

limma_de2=reshape2::dcast(limma_de,pos+SYMBOL+annotation1~group)
table(limma_de2$Hypomethylation>0)
table(limma_de2$Hypomethylation>0,limma_de2$Hypermethylation>0)
table(limma_de2$Hypomethylation>10,limma_de2$Hypermethylation>10)
limma_de2$group=limma_de[match(limma_de2$pos,limma_de$pos),"group"]
limma_de2=limma_de2[limma_de2$Hypermethylation>10 | limma_de2$Hypomethylation>10,]

hgnc_complete_set <- read.delim("~/database/hgnc_complete_set.txt")
table(limma_de2$SYMBOL %in% hgnc_complete_set$symbol)

limma_de2$type=hgnc_complete_set[match(limma_de2$SYMBOL,hgnc_complete_set$symbol),"locus_group"]
limma_de2$type=ifelse( !is.na(limma_de2$type),limma_de2$type,gtf22[match(limma_de2$SYMBOL,gtf22$gene_name),"gene_type"])
g=unique(limma_de2[is.na(limma_de2$type),"SYMBOL"])

unique(limma_de2$type)
limma_de2$type=ifelse(limma_de2$type %in% c("protein-coding gene","protein_coding"),"protein-coding gene",
                      ifelse(limma_de2$type %in% c("lncRNA","lincRNA","non-coding RNA"),"non-coding RNA",
                             ifelse(limma_de2$type %in% c("pseudogene","unprocessed_pseudogene"),"pseudogene",limma_de2$type)))
write.table(limma_de2,"limma_de_top1000.tsv",quote = F,sep = "\t",row.names = F)

table(limma_de2$type)
table(limma_de2$group)

table(limma_de2$type,limma_de2$group)
table(limma_de2$annotation1,limma_de2$group)
p_value1=chisq.test(table(limma_de2$type,limma_de2$group))$p.value
p_value2=chisq.test(table(limma_de2$annotation1,limma_de2$group))$p.value
a= limma_de2[,c("type","group")] %>% dplyr::group_by(group,type) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
b= limma_de2[,c("group","annotation1")] %>% dplyr::group_by(group,annotation1) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
a$n_lab=ifelse(a$Freq<1,"",a$n)
a$Freq_lab=ifelse(a$Freq<1,"",paste(round(a$Freq,2),"%"))
b$n_lab=ifelse(b$Freq<1,"",b$n)
b$Freq_lab=ifelse(b$Freq<1,"",paste(round(b$Freq,2),"%"))

sfig_2a=ggplot(a[!is.na(a$type),],aes(x=group,y=Freq,fill=type))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))
sfig_2b=ggplot(b,aes(x=group,y=Freq,fill=annotation1))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))


## lasso               ##########################
rm(list=ls())
gc()

load("lasso_gene21.RData")
plot(cv_fit) # fig_2a

## rf important        #########################################
load("medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21.RData")
all_imp=data.frame()
for (i in 1:100) {
    all_imp=rbind(all_imp,data.frame(pos=rownames(limma_modlist[[i]]$rf_mod$finalModel$importance),importance=limma_modlist[[i]]$rf_mod$finalModel$importance))
}

a=do.call(rbind,lapply(unique(all_imp$pos), function(x){data.frame(pos=x,med=median(all_imp[all_imp$pos %in% x,"MeanDecreaseGini"]))}))
a=a[order(a$med,decreasing = T),]
all_imp$pos=factor(all_imp$pos,levels = a$pos)

load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
imp_de_gene2=imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,]
imp_de_gene2$group=imp_de_gene[match(imp_de_gene2$pos,imp_de_gene$pos),"group"]

hgnc_complete_set <- read.delim("~/database/hgnc_complete_set.txt")
imp_de_gene2$type=hgnc_complete_set[match(imp_de_gene2$SYMBOL,hgnc_complete_set$symbol),"locus_group"]
all_imp=merge(all_imp,imp_de_gene2,by="pos")

fig_2b=ggplot(all_imp, aes(pos, MeanDecreaseGini)) +
    geom_boxplot(aes(fill=group))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    scale_fill_brewer(palette = "Set1")+
    labs(title = "")+xlab(label = "")+ylab("MeanDecreaseGini")+theme(axis.text.x = element_blank())


## AUC Cohort          ##########################################
rm(list=ls())
load("medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21.RData")
load("medip_stad_peak_samp250_subsamp.RData")

all_test_pre=data.frame()
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    all_test_pre=rbind(all_test_pre,limma_pre[limma_pre$sample %in% test_samp,])
    all_valid_pre=rbind(all_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}

g1=pROC::ggroc(list(`Discovery set`=roc(as.factor(all_test_pre$group), as.numeric(all_test_pre$pre)),
                    `Validation set`=roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre))))
auc1=roc(as.factor(all_test_pre$group), as.numeric(all_test_pre$pre))
auc2=roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre))
fig_2c=g1+ggsci::scale_color_lancet()+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    theme_classic() + theme(legend.position = "")+
    annotate("text", y = 0.5, x=0.5, label = paste0("Discovery set\n AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )",
                                                    "\n\nValidation set\n AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )"),colour ="#00468BFF")

## phen risk roc test  ########################
load("medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21.RData")
load("medip_stad_peak_samp250_subsamp.RData")
all_test_pre=data.frame()
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    all_test_pre=rbind(all_test_pre,limma_pre[limma_pre$sample %in% test_samp,])
    all_valid_pre=rbind(all_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}
all_pre=rbind(all_test_pre,all_valid_pre)

all_pre=do.call(rbind,lapply(unique(all_pre$sample),function(x){
    data.frame(sample=x,risk=median(all_pre[all_pre$sample %in% x,"pre"]))
}))

load("medip_stad_peak_samp250_phen.RData")
table(phen$Stage)
phen=merge(phen,all_pre,by="sample")

for (i in c("III","I","IV","II","0")) {
    dat=phen[phen$Stage %in% c("normal",i),]
    assign(paste0("Stage",i,"_auc"),roc(as.factor(dat$Group), as.numeric(dat$risk)))
    assign(paste0("Stage",i,"_auc1"),roc(as.factor(dat$Group), as.numeric(dat$risk))$auc[[1]])
}
g1=pROC::ggroc(list(`Stage 0`=Stage0_auc,`Stage I`=StageI_auc,`Stage II`=StageII_auc,
                    `Stage III`=StageIII_auc,`Stage IV`=StageIV_auc))

fig_2d=g1+ggsci::scale_color_lancet() + labs(x = "1 - Specificity", y = "Sensitivity") + 
    theme_classic() + theme(legend.title = element_blank())+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    annotate("text", y = 0.5, x=0.5, label = paste0("Stage 0 AUC = ",round(Stage0_auc1,4)," \nStage I AUC = ",round(StageI_auc1,4)," \nStage II AUC = ",round(StageII_auc1,4)," \nStage III AUC = ",round(StageIII_auc1,4)," \nStage IV AUC = ",round(StageIV_auc1,4)),colour ="#00468B99")



g1=pROC::ggroc(list(Age=roc(as.factor(phen$Group), as.numeric(phen$Age)),
                    Gender=roc(as.factor(phen$Group), as.numeric(phen$Gender)),
                    Concentration=roc(as.factor(phen$Group), as.numeric(phen$Concentration)),
                    Risk=roc(as.factor(phen$Group), as.numeric(phen$risk))))
auc1=roc(as.factor(phen$Group), as.numeric(phen$Age))$auc[[1]]
auc2=roc(as.factor(phen$Group), as.numeric(phen$Gender))$auc[[1]]
auc3=roc(as.factor(phen$Group), as.numeric(phen$Concentration))$auc[[1]]
auc4=roc(as.factor(phen$Group), as.numeric(phen$risk))$auc[[1]]


fig_2e=g1+ggsci::scale_color_lancet() + labs(x = "1 - Specificity", y = "Sensitivity") + 
    theme_classic() + theme(legend.title = element_blank())+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    annotate("text", y = 0.5, x=0.5, label = paste0("                 Age AUC = ",round(auc1,4)," \n           Gender AUC = ",round(auc2,4)," \nConcentration AUC = ",round(auc3,4)," \n               Risk AUC = ",round(auc4,4)),colour ="#00468B99")



median(phen$Age,na.rm = T)
phen$Age=ifelse(phen$Age<56,"<56",ifelse(phen$Age>=56,">=56",NA))
median(phen$Concentration,na.rm = T)
phen$Concentration=ifelse(phen$Concentration<7.25,"<7.25",ifelse(phen$Concentration>=7.25,">=7.25",NA))

dat=reshape2::melt(phen[,c(1:5,8)],c("sample","risk","Group"))
dat$sub_group=paste0(dat$variable," ",dat$value)
dat$risk=1-dat$risk
dat=dat[!is.na(dat$value),]

dat2=reshape2::melt(phen[,c(1,5,6,8)],c("sample","risk","Group"))
dat2=dat2[!is.na(dat2$value),]
dat2$sub_group=ifelse(!dat2$value %in% "normal",paste0(dat2$variable," ",dat2$value),dat2$value)
dat2$risk=1-dat2$risk

dat3=reshape2::melt(phen[,c(1,5,6,8)],c("sample","risk","Group"))
dat3=dat3[!is.na(dat3$value),]
dat3=dat3[!dat3$value %in% "normal",]
dat3$sub_group=ifelse(!dat3$value %in% "normal",paste0(dat3$variable," ",dat3$value),dat3$value)
dat3$risk=1-dat3$risk


sfig_3a=ggplot(dat, aes(sub_group, risk,fill=sub_group)) + 
    geom_boxplot()+geom_jitter(width = 0.1, alpha = 0.5, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    stat_compare_means(aes(group=sub_group),label.y = 1.1)+
    facet_grid(.~Group)+theme_bw()+ylab("Predicted risk score")+
    #theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),panel.grid = element_blank())
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),legend.title = element_blank())

sfig_3b=ggplot(dat3, aes(sub_group, risk,fill=sub_group)) + 
    geom_boxplot()+geom_jitter(width = 0.1, alpha = 0.5, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    stat_compare_means(aes(group=sub_group),label.y = 1.1)+
    theme_bw()+ylab("Predicted risk score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),legend.title = element_blank())


## with cfdna          ##################################################
load("medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21_with_cfdna.RData")
limma_modlist_cfdna=limma_modlist

load("medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21.RData")
load("medip_stad_peak_samp250_subsamp.RData")

all_test_pre=data.frame()
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    all_test_pre=rbind(all_test_pre,limma_pre[limma_pre$sample %in% test_samp,])
    all_valid_pre=rbind(all_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}
all_pre=rbind(all_test_pre,all_valid_pre)

cfdna_test_pre=data.frame()
cfdna_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist_cfdna[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    cfdna_test_pre=rbind(cfdna_test_pre,limma_pre[limma_pre$sample %in% test_samp,])
    cfdna_valid_pre=rbind(cfdna_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}
cfdna_pre=rbind(cfdna_test_pre,cfdna_valid_pre)

library(pROC)
library(rms)
rcorrcens(cfdna_pre$group~cfdna_pre$pre)
ci(roc(as.factor(cfdna_pre$group), as.numeric(cfdna_pre$pre)))
coords(roc(as.factor(cfdna_pre$group), as.numeric(cfdna_pre$pre)), "best", ret = "all", transpose = FALSE)
roc(as.factor(cfdna_pre$group), as.numeric(cfdna_pre$pre))


library(rms)
g1=pROC::ggroc(list(`RF Model`=roc(as.factor(all_pre$group), as.numeric(all_pre$pre)),
                    `Combine Model`=roc(as.factor(cfdna_pre$group), as.numeric(cfdna_pre$pre))
))
auc1=roc(as.factor(all_pre$group), as.numeric(all_pre$pre))
auc2=roc(as.factor(cfdna_pre$group), as.numeric(cfdna_pre$pre))

fig_2f=g1+ggsci::scale_color_lancet()+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    theme_classic() + theme(legend.title = element_blank())+
    annotate("text", y = 0.75, x=0.5, label = paste0("p=",format(roc.test(auc1, auc2)$p.value,scientific=T,digits=4)),colour ="#00468BFF")+
    annotate("text", y = 0.6, x=0.5, label = paste0("RF Model AUC = ",round(auc1$auc[[1]],4),"\n95% CI ( ",round(ci(auc1)[1],4)," - ",round(ci(auc1)[3],4)," )"),colour ="#00468BFF")+
    annotate("text", y = 0.45, x=0.5, label = paste0("\nCombine Model AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )"),colour ="#ED0000FF")

