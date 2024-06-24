rm(list=ls())
options(stringsAsFactors = F)

load("medip_stad_peak_samp250_phen.RData")
all_coverage <- read.table("all_coverage.txt", quote="\"", comment.char="")
all_coverage$V1=basename(all_coverage$V1)
all_coverage$V1=gsub("-","_",all_coverage$V1)
all_coverage$V1=gsub("_samtools_rmdup_sorted.bam","",all_coverage$V1)
table(phen$sample %in% all_coverage$V1)
all_coverage=all_coverage[all_coverage$V1 %in% phen$sample,]
save(all_coverage,file = "peak_all_coverage.RData")


load("peak_all_coverage.RData")
load("peak_raw_exp.RData")
load("medip_stad_peak_samp250_phen.RData")

train_test=phen[phen$group == "train_test",]
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
save(samp_select,file = "samp_select.RData")


load("samp_select.RData")
train_test=phen[phen$group %in% "train_test" & phen$sample %in% samp_select,]
train_test_samp=train_test$sample

train_test=peak_exp[,train_test_samp]
keep <- rowSums(train_test > 0) >= ncol(train_test)*0.2
table(keep)
save(keep,file = "peak_keep.RData")

train_test=train_test[keep,]
save(train_test,file="peak_train_test_qc_data.RData")

f75 <- rep_len(1,ncol(train_test))
for (j in seq_len(ncol(train_test))) f75[j] <- quantile(train_test[,j], probs=0.75)
lib.size=colSums(train_test)
f75=f75 / lib.size
refColumn <- colnames(train_test)[which.min(abs(f75-mean(f75)))]
save(refColumn,file="peak_refColumn.RData")

valid=phen[phen$group %in% "valid" & phen$sample %in% samp_select,]
valid_samp=valid$sample
valid=peak_exp[keep,c(valid_samp,refColumn)]
save(valid,file="peak_valid_qc_data.RData")




