## medip model com     #####################################################
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
save(limma_modlist,file="medip_stad_peak_samp250_norm_voom_10_02_limma_modlist_PSM_gene3_cv_com_lasso_de21_with_cfdna2.RData")


