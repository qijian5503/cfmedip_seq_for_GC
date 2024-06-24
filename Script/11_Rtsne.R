## Rtsne               ###########################################################
load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
table(imp_de_gene2$Hypermethylation>10)
table(imp_de_gene2$Hypomethylation>10)
de_pos=as.character(imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,"pos"])

load("traintest_norm.RData")
load("valid_norm.RData")
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

p1=ggplot(dat_plt,aes(x=dat_plt[,1],y=dat_plt[,2],color=Group))+geom_point()+
    scale_color_manual(values=c("#0571B0","#CA0020"))+
    theme(legend.title = element_text(size = 13),legend.text = element_text(size = 13),axis.title.y= element_text(size=15),axis.title.x= element_text(size=15),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+
    labs(x = "",y = "",title ="" )
p1

