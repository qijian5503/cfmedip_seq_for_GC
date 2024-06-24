load("limma_modlist_PSM")
load("medip_stad_peak_samp250_subsamp.RData")

library(pROC)
all_test_pre=data.frame()
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    test_pre=limma_pre[limma_pre$sample %in% test_samp,]
    valid_pre=limma_pre[!limma_pre$sample %in% test_samp,]
    
    test_auc=roc(as.factor(test_pre$group), as.numeric(test_pre$pre))$auc[1]
    valid_auc=roc(as.factor(valid_pre$group), as.numeric(valid_pre$pre))$auc[1]
    test_pre=coords(roc(as.factor(test_pre$group), as.numeric(test_pre$pre)), "best", ret = "all", transpose = FALSE)
    valid_pre=coords(roc(as.factor(valid_pre$group), as.numeric(valid_pre$pre)), "best", ret = "all", transpose = FALSE)
    
    all_test_pre=rbind(all_test_pre,data.frame(auc=test_auc,threshold=test_pre[[1]],specificity=test_pre[[2]],sensitivity=test_pre[[3]],accuracy=test_pre[[4]]))
    all_valid_pre=rbind(all_valid_pre,data.frame(auc=valid_auc,threshold=valid_pre[[1]],specificity=valid_pre[[2]],sensitivity=valid_pre[[3]],accuracy=valid_pre[[4]]))
}

write.table(all_test_pre,"test_set_cutoff_auc.tsv",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(all_valid_pre,"valid_set_cutoff_auc.tsv",quote = F,row.names = F,col.names = F,sep = "\t")


pdf("rf_tree_plot.pdf",20,10)
for (i in 1:100) {
    for (k in 1:100) {
        p=reprtree::plot.getTree(limma_modlist[[i]]$rf_mod,k = k,caret = T)
        print(p)
    }
}
dev.off()
