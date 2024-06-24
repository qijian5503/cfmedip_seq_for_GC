## medip model voom    #####################################################
rm(list=ls())
load("limma_PSM_de.RData")
load("peakAnno.RData")
limma_de=list()
for(i in 1:100){
    res=limma_de_all[[i]]
    limma_de[[i]]=res[abs(res$logFC)> 3 & res$adj.P.Val<0.0001,]
    message(nrow(limma_de[[i]]))
}

imp_de_gene=data.frame()
for(i in 1:100){
    de=limma_de[[i]]
    de$pos=rownames(de)
    imp_de_gene=rbind(imp_de_gene,de[,c(1,5,7)])
}
rownames(imp_de_gene)=NULL

imp_de_gene$group=ifelse(imp_de_gene$logFC>0,"Hypermethylation","Hypomethylation")
imp_de_gene=merge(imp_de_gene,peak_ana[,c(1,7,8)],by="pos",all.x=T)
#save(imp_de_gene,file = "imp_de_gene.RData")


load("medip_stad_peak_samp250_subsamp.RData")
load("medip_stad_peak_samp250_phen.RData")
load("samp_select.RData")
load("traintest_norm.RData")
load("valid_norm.RData")
valid_sample=colnames(norm_v)

load("imp_de_gene.RData")
de_pos=as.data.frame(table(imp_de_gene$pos))
table(de_pos$Freq>10)
de_pos=as.character(de_pos[de_pos$Freq>10 ,1])

library(glmnet)
de_exp=as.data.frame(t(norm[de_pos,]))
de_exp$Group=phen[rownames(de_exp),"Group"]

cv_fit = cv.glmnet(as.matrix(de_exp[,de_pos]), de_exp$Group, family = "binomial",type.measure = "class")
save(cv_fit,file = "lasso_gene21.RData")
coef.min = coef(cv_fit, s = "lambda.min")
de_gene=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
plot(cv_fit)

library(glmnet)
library(randomForest)
ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% samp_select]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% samp_select]
    de_exp=as.data.frame(t(norm[de_gene,]))
    de_exp$group=factor(phen[rownames(de_exp),"Group"],levels = c("cancer","normal"))
    train_data=de_exp[train_samp,c(de_gene,"group")]
    test_data=de_exp[test_samp,c(de_gene,"group")]
    
    rf_mod = train(group ~ .,data = train_data,method = "rf",trControl = ctrl)
    rf_pre=predict(rf_mod,test_data,type="prob")%>%data.frame
    rf_auc=unlist(slot(performance(prediction(rf_pre$normal,test_data$group),'auc'),"y.values"))
    valid_data=as.data.frame(t(norm_v[de_gene,]))
    valid_pre=predict(rf_mod,valid_data,type="prob")%>%data.frame
    valid_pre$group=as.factor(phen[rownames(valid_pre),"Group"])
    valid_auc=unlist(slot(performance(prediction(valid_pre$normal,valid_pre$group),'auc'),"y.values"))
    
    limma_pre=rbind(data.frame(sample=rownames(test_data),group=test_data$group,pre=rf_pre$normal),data.frame(sample=rownames(valid_pre),group=valid_pre$group,pre=valid_pre$normal))
    limma_modlist[[i]]= list(limma_pre=limma_pre,rf_mod=rf_mod,limma_de=de,valid_auc=valid_auc,rf_auc=rf_auc)
    message(round(rf_auc,3),"   ",round(valid_auc,3))
}
save(limma_modlist,file="limma_modlist_PSM_rf.RData") 

