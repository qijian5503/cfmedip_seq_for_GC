load("peak_train_test_qc_data.RData")
d=DGEList(counts = train_test) %>% calcNormFactors(. ,method = "TMM" ,refColumn = refColumn) %>% voom()
norm=as.data.frame(d$E)
norm[1:3,1:3]
save(norm,file = "traintest_norm.RData")

load("peak_valid_qc_data.RData")
d=DGEList(counts = valid) %>% calcNormFactors(. ,method = "TMM" ,refColumn = refColumn) %>% voom()
norm_v=as.data.frame(d$E)
table(norm_v[,refColumn]==norm[,refColumn])
norm_v[,refColumn]=NULL
save(norm_v,file = "valid_norm.RData")

ref_exp=data.frame(V1=rownames(valid),V2=valid[,refColumn])
save(ref_exp,file = "ref_exp.RData")



