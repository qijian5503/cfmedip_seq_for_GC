## createDataPartition  ############################
load("medip_stad_peak_samp250_phen.RData")
table(phen$group)
train_test=phen[phen$group %in% "train_test",]
train_test_samp=train_test$sample
subsample=createDataPartition(train_test$Group, p = 0.8,times = 100)
for (i in 1:100) {
    train_samp=train_test_samp[subsample[[i]]]
    test_samp=train_test_samp[!train_test_samp %in% train_samp]
    subsample[[i]]=list(train_samp=train_samp,test_samp=test_samp)
}
save(subsample,file="medip_stad_peak_samp250_subsamp.RData")

