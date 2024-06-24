rm(list=ls())
gc()
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

library(rms)
library(pROC)
coords(roc(as.factor(all_pre$group), as.numeric(all_pre$pre)), "best", ret = "all", transpose = FALSE)
# threshold specificity sensitivity  accuracy   tn   tp  fn  fp
# 0.465   0.9218314   0.9309357 0.9251908 4953 2925 217 420

all_pre$pre_group=ifelse(all_pre$pre>0.465,"normal","cancer")
a=as.data.frame(table(all_pre$sample,all_pre$pre_group))
a=reshape2::dcast(a,Var1~Var2,value.var = "Freq")
a$group2=ifelse(a$normal>a$cancer,"normal","cancer")
a$group=all_pre[match(a$Var1,all_pre$sample),"group"]
table(a$group,a$group2)
a$Classification=ifelse(a$group2 == a$group,"Correct","Incorrect")
table(a$Classification)

load("medip_stad_peak_samp250_phen.RData")
load("peak_tab.RData")
peak_tab$Classification=a[match(peak_tab$Sample,a$Var1),"Classification"]
peak_tab$`Age(years)`=phen[match(peak_tab$Sample,phen$Sample),"Age"]

library(gtsummary)
colnames(peak_tab)
tbl_summary(peak_tab[,c(2:6,8:11,32,33)],by = Classification,
            statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
            digits = all_continuous() ~ 2,#确定小数点数
            missing_text = "(Missing)") %>% add_p() %>%  add_overall() %>% bold_labels()


