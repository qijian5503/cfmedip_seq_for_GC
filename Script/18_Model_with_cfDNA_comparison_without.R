load("limma_modlist_PSM_lasso_with_cfdna.RData")
limma_modlist_cfdna=limma_modlist

load("limma_modlist_PSM_rf.RData")
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


library(rms)
g1=pROC::ggroc(list(`RF Model`=roc(as.factor(all_pre$group), as.numeric(all_pre$pre)),
                    `Combine Model`=roc(as.factor(cfdna_pre$group), as.numeric(cfdna_pre$pre))
))
auc1=roc(as.factor(all_pre$group), as.numeric(all_pre$pre))
auc2=roc(as.factor(cfdna_pre$group), as.numeric(cfdna_pre$pre))
round(ci(auc1)[1],4)
round(ci(auc2)[3],4)

p1=g1+ggsci::scale_color_lancet()+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    theme_classic() + theme(legend.title = element_blank())+
    annotate("text", y = 0.75, x=0.5, label = paste0("p=",format(roc.test(auc1, auc2)$p.value,scientific=T,digits=4)),colour ="#00468BFF")+
    annotate("text", y = 0.6, x=0.5, label = paste0("RF Model AUC = ",round(auc1$auc[[1]],4),"\n95% CI ( ",round(ci(auc1)[1],4)," - ",round(ci(auc1)[3],4)," )"),colour ="#00468BFF")+
    annotate("text", y = 0.45, x=0.5, label = paste0("\nCombine Model AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )"),colour ="#ED0000FF")
p1


