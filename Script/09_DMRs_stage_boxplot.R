## single pos boxplot  ########################
load("lasso_gene21.RData")
coef.min = coef(cv_fit, s = "lambda.min")
de_gene=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因

load("traintest_norm.RData")
load("valid_norm.RData")
dat=cbind(norm,norm_v)
dat=as.data.frame(t(dat[de_gene,]))

load("medip_stad_peak_samp250_phen.RData")
dat=cbind(dat,phen[rownames(dat),])

library(ggpubr)
pdf("fig_s3.pdf")
for (i in de_gene) {
    dat$pos=dat[,i]
    p=ggplot(dat, aes(Group, pos)) +
        geom_boxplot(aes(fill=Group))+ geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
        stat_compare_means(aes(group=Group),hide.ns = T)+scale_fill_brewer(palette = "Set1")+
        labs(title = "")+xlab(label = "")+ylab(paste0(i," Methylation"))
    print(p)
}
dev.off()

