## TCGA rnaseq expression       ##################################
rm(list=ls())
load("lasso_gene21.RData")
coef.min = coef(cv_fit, s = "lambda.min")
de_gene=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
load("imp_de_gene.RData")
imp_de_gene=imp_de_gene[imp_de_gene$pos %in% de_gene,]
imp_de_gene=unique(imp_de_gene[,c(1,5,6,7)])


load("tcga_gtex_phenotype.RData")
load("tcga_all_pheno.RData")
load("gencode_v23_annotation_gtf.RData")
dat=data.table::fread("TcgaTargetGtex_rsem_gene_tpm.gz")
dat=as.data.frame(dat)
imp_de_gene$SYMBOL[!imp_de_gene$SYMBOL %in% gtf$gene_name]

gtf=gtf[gtf$gene_name %in% imp_de_gene$SYMBOL,]
dat=merge(gtf[,-2],dat,by.x="gene_id",by.y="sample")
dat=limma::avereps(dat[,3:ncol(dat)],ID=dat$gene_name)
dat=as.data.frame(t(dat))

table(tcga_gtex_phenotype$sample %in% rownames(dat))
samp=tcga_gtex_phenotype$sample
gen_exp=as.data.frame(dat[samp,])
anyNA(gen_exp)
gen_exp$sample=rownames(gen_exp)
gen_exp=merge(gen_exp,tcga_gtex_phenotype[,c(2,8,9)],by="sample")
table(gen_exp$tcga_type)
gen_exp$Tumor_type=factor(gen_exp$Tumor_type,levels = c("Tumor","Normal"))

load("tcga_gtex_surv.RData")
table(gen_exp$sample %in% rownames(surv))
gen_exp$OS.time=surv[match(gen_exp$sample,rownames(surv)),"OS.time"]
gen_exp$OS=surv[match(gen_exp$sample,rownames(surv)),"OS"]
save(gen_exp,file = "tcga_gen_exp.RData")


load("tcga_gen_exp.RData")
min(gen_exp[,2:15])
gen_exp[,2:15]=2**(gen_exp[,2:15])
gen_exp[,2:15]=gen_exp[,2:15]-min(gen_exp[,2:15])
gen_exp[,2:15]=log2(gen_exp[,2:15]+1)

library(pROC)
can=unique(gen_exp$tcga_type)[!unique(gen_exp$tcga_type) %in% c("UVM","MESO","SARC","PCPG")]

all_res=data.frame()
for (g in colnames(gen_exp)[2:15]) {
    dat1=gen_exp[,c("Tumor_type","tcga_type",g)]
    de=do.call(rbind,lapply(can, function(x){
        data.frame(cancer=x,
                   lfc=log(median(as.numeric(dat1[dat1$tcga_type %in% x & dat1$Tumor_type %in% "Tumor",g]),na.rm = T)+0.00001)-log(median(as.numeric(dat1[dat1$tcga_type %in% x & dat1$Tumor_type %in% "Normal",g]),na.rm = T)+0.00001),
                   auc=roc(dat1[dat1$tcga_type %in% x,"Tumor_type"], as.numeric(dat1[dat1$tcga_type %in% x,g]))$auc,
                   de_p=wilcox.test(dat1[dat1$tcga_type %in% x,g]~dat1[dat1$tcga_type %in% x,"Tumor_type"])$p.value)
    }))
    all_res=rbind(all_res,data.frame(gene=g,de))
}
all_res$de_p=as.numeric(all_res$de_p)
all_res$auc=as.numeric(all_res$auc)
save(all_res,file = "tcga_exp_all_res.RData")


library(survival)
library(survminer)
load("tcga_gen_exp.RData")
surv=gen_exp[!is.na(gen_exp$OS),]
gen_surv_result=data.frame()
for (g in colnames(gen_exp)[2:15]) {
    for (i in unique(surv$tcga_type)) {
        dat=surv[surv$tcga_type==i,]
        result = tryCatch({
            dat$exp=dat[,g]
            dat$group=ifelse(dat[,g] > median(dat[,g],na.rm = T),"High","Low")
            km_p=surv_pvalue(survfit(Surv(OS.time, OS)~group,data = dat))$pval
            cox_res=summary(coxph(Surv(OS.time, OS)~exp,data = dat))
            data.frame(gene=g,tcga_type=i,samp_num=nrow(dat),km_p=km_p,cox_p=cox_res[["coefficients"]][5],cox_hr=cox_res[["conf.int"]][1],low_CI=cox_res[["conf.int"]][3],up_CI=cox_res[["conf.int"]][4])
        }, error = function(e){data.frame(gene=g,tcga_type=i,samp_num=nrow(dat),km_p="",cox_p="",cox_hr="",low_CI="",up_CI="")})
        gen_surv_result=rbind(gen_surv_result,result)
    }
}
save(gen_surv_result,file = "tcga_exp_surv_result.RData")



## TCGA methylation             ##################################
load("tcga_gtex_phenotype.RData")
can_nam=unique(tcga_gtex_phenotype$tcga_type)
library(data.table)
library(ChAMP)
library(impute)
for (i in can_nam) {
    if(file.exists(paste0("/TCGA-",i,".methylation450.tsv"))){
        methy=fread(paste0("/TCGA-",i,".methylation450.tsv"),stringsAsFactors = F,data.table = F)
    }else{methy=fread(paste0("/TCGA-",i,".methylation450.tsv.gz"),stringsAsFactors = F,data.table = F)}
    rownames(methy)=methy$`Composite Element REF`
    methy$`Composite Element REF`=NULL
    methy=methy[rowSums(is.na(methy))<ncol(methy)*0.5,]
    methy=impute.knn(as.matrix(methy))
    condition=data.frame(row.names = colnames(methy$data), sample=colnames(methy$data),group=substr(colnames(methy$data),14,16))
    methy=champ.filter(beta = methy$data ,pd = condition)
    methy=champ.norm(beta=methy$beta,arraytype="450K",cores=8)
    save(methy,file=paste0(i,"_methy_450k_norm.RData"))
}


load("lasso_gene21.RData")
coef.min = coef(cv_fit, s = "lambda.min")
de_gene=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
load("imp_de_gene.RData")
imp_de_gene=imp_de_gene[imp_de_gene$pos %in% de_gene,]
imp_de_gene=unique(imp_de_gene[,c(1,5,6,7)])

gpl=read.delim("Annotation_file/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy", header=FALSE, comment.char="#")
table(imp_de_gene$SYMBOL %in% gpl$V2)
id=unique(as.character(gpl[gpl$V2 %in% imp_de_gene$SYMBOL,1]))

load(paste0(can_nam[1],"_methy_450k_norm.RData"))
all_methy_data=as.data.frame(methy[rownames(methy) %in% id,])
all_methy_data$id=rownames(all_methy_data)
for (i in can_nam[-1]) {
    load(paste0(i,"_methy_450k_norm.RData"))
    dat=as.data.frame(methy[rownames(methy) %in% id,])
    dat$id=rownames(dat)
    all_methy_data=merge(all_methy_data,dat,by="id",all=T)
}
rownames(all_methy_data)=all_methy_data$id
all_methy_data$id=NULL
all_methy_data=as.data.frame(t(all_methy_data))
all_methy_data$sample=substr(rownames(all_methy_data),1,15)
dat=merge(all_methy_data,tcga_gtex_phenotype[,c(2,8,9)],by="sample")

load("tcga_gtex_surv.RData")
table(dat$sample %in% rownames(surv))
dat$OS.time=surv[match(dat$sample,rownames(surv)),"OS.time"]
dat$OS=surv[match(dat$sample,rownames(surv)),"OS"]
save(dat,file = "methy_data.RData")



load("methy_data.RData")
min(dat[,2:30])
a=table(dat$tcga_type,dat$Tumor_type)%>% as.data.frame()
a=a[a$Var2 %in% "Normal" & a$Freq>5,]
can_nam=as.character(a$Var1)
dat=dat[dat$tcga_type %in% can_nam,]

gpl=read.delim("E:/database/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy", header=FALSE, comment.char="#")
gpl=gpl[gpl$V1 %in% colnames(dat),]
table(gpl$V2)

all_dat=cbind(dat[,c(1,31:34)],data.frame(ANKRD30BL=as.numeric(rowMeans(dat[,gpl[gpl$V2 %in% "ANKRD30BL","V1"]])), 
                                          CCDC144B=as.numeric(dat[,"cg10441789"]), 
                                          PKP4=as.numeric(rowMeans(dat[,gpl[gpl$V2 %in% "PKP4","V1"]])), 
                                          RRN3P1=as.numeric(rowMeans(dat[,gpl[gpl$V2 %in% "RRN3P1","V1"]])), 
                                          ZNF74=as.numeric(rowMeans(dat[,gpl[gpl$V2 %in% "ZNF74","V1"]]))))

all_res=data.frame()
for (g in colnames(all_dat)[6:10]) {
    dat1=all_dat[,c("Tumor_type","tcga_type",g)]
    de=do.call(rbind,lapply(can_nam, function(x){
        data.frame(cancer=x,
                   lfc=log(median(as.numeric(dat1[dat1$tcga_type %in% x & dat1$Tumor_type %in% "Tumor",g]),na.rm = T))-log(median(as.numeric(dat1[dat1$tcga_type %in% x & dat1$Tumor_type %in% "Normal",g]),na.rm = T)),
                   auc=roc(dat1[dat1$tcga_type %in% x,"Tumor_type"], as.numeric(dat1[dat1$tcga_type %in% x,g]))$auc,
                   de_p=wilcox.test(dat1[dat1$tcga_type %in% x,g]~dat1[dat1$tcga_type %in% x,"Tumor_type"])$p.value)
    }))
    all_res=rbind(all_res,data.frame(id=g,de))
}
all_res$de_p=as.numeric(all_res$de_p)
all_res$auc=as.numeric(all_res$auc)
save(all_res,file = "tcga_methy_all_res.RData")

library(survival)
library(survminer)
load("methy_data.RData")
surv=all_dat[!is.na(all_dat$OS.time),]
methy_surv_result=data.frame()
for (g in colnames(all_dat)[6:10]) {
    for (i in unique(surv$tcga_type)) {
        dat=surv[surv$tcga_type==i,]
        result = tryCatch({
            dat$exp=dat[,g]
            dat$group=ifelse(dat[,g] > median(dat[,g],na.rm = T),"High","Low")
            km_p=surv_pvalue(survfit(Surv(OS.time, OS)~group,data = dat))$pval
            cox_res=summary(coxph(Surv(OS.time, OS)~exp,data = dat))
            data.frame(gene=g,tcga_type=i,samp_num=nrow(dat),km_p=km_p,cox_p=cox_res[["coefficients"]][5],cox_hr=cox_res[["conf.int"]][1],low_CI=cox_res[["conf.int"]][3],up_CI=cox_res[["conf.int"]][4])
        }, error = function(e){data.frame(gene=g,tcga_type=i,samp_num=nrow(dat),km_p="",cox_p="",cox_hr="",low_CI="",up_CI="")})
        methy_surv_result=rbind(methy_surv_result,result)
    }
}
save(methy_surv_result,file = "tcga_methy_surv_result.RData")

## tcga rnaseq expression plot  #####################################
rm(list=ls())
load("tcga_exp_all_res.RData")
a=all_res[all_res$auc>0.8 | all_res$auc<0.2,]
b=all_res[all_res$de_p<0.05,]
table(a$gene,a$cancer)
table(b$gene,b$cancer)

colnames(all_res)
a=reshape2::dcast(all_res,cancer~gene,value.var = "auc")
summary(all_res$auc)
rownames(a)=a$cancer
a$cancer=NULL
a=a[,colSums(a>0.8)>1]
library(ComplexHeatmap)
library(circlize)
col_fun = circlize::colorRamp2(c(0,0.5, 1), c("#0571B0","white", "#CA0020"))
star_p=function(x){if(x>=0.8){format(x,scientific=F,digits=2)}else{" "}}
p1 = Heatmap(a, name = "AUC",column_title = "",col = col_fun,row_names_side = "left",column_names_side = "bottom", rect_gp = gpar(type = "none"),border = TRUE, 
             cell_fun = function(j, i, x, y, width, height, fill) {
                 if(a[i, j] != 0){
                     grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = col_fun(a[i, j]), fill =  col_fun(a[i, j])))
                     grid.text(star_p(a[i, j]), x, y, gp = gpar(fontsize = 8))
                 }else{grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill =  "white"))}}, cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = T, show_column_names = T)

colnames(all_res)
all_res=all_res[!is.na(all_res$auc),]
all_res$change = as.factor(ifelse(all_res$lfc> 0 ,'Up','Down'))
all_res=all_res[!is.na(all_res$change),]
all_res$`-log10(P value)`=-log10(all_res$de_p)
p2=ggplot(na.omit(all_res[all_res$de_p<0.05,]), aes_string(x="gene", y="cancer", size="`-log10(P value)`", color="change")) +
    geom_point(aes())+ xlab("")+ylab("")+ 
    scale_color_manual(values = c("#CA0020","#0571B0"))+
    scale_size()+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
p2

load("tcga_exp_surv_result.RData")
gen_surv_result$cox_p=as.numeric(gen_surv_result$cox_p)
gen_surv_result$cox_hr=as.numeric(gen_surv_result$cox_hr)
gen_surv_result=gen_surv_result[!is.na(gen_surv_result$cox_p),]
dat=gen_surv_result[gen_surv_result$cox_p<0.05,]
dat$`-log10(COX P value)`=-log10(dat$cox_p)
dat$`-log10(COX HR)`=log10(dat$cox_hr)
p3=ggplot(dat, aes_string(y="tcga_type", x="gene", color="`-log10(COX HR)`", size="`-log10(COX P value)`")) +
    geom_point(aes())+ xlab("")+ylab("")+ scale_color_gradient(high="#CA0020",low="#0571B0")+
    scale_size()+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


## tcga methylation plot        #####################################
load("tcga_methy_all_res.RData")
table(all_res$auc>0.8 | all_res$auc<0.2)
a=all_res[all_res$auc>0.8 | all_res$auc<0.2,]
table(a$id,a$cancer)

b=all_res[all_res$de_p<0.05,]
table(b$id,b$cancer)
table(b$id)
table(b$cancer)

colnames(all_res)
all_res$change = as.factor(ifelse(all_res$lfc> 0 ,'Hypermethylation','Hypomethylation'))
all_res=all_res[!is.na(all_res$change),]
all_res$LogFC=all_res$lfc
all_res$`-log10(P value)`=-log(all_res$de_p)

p4=ggplot(all_res[all_res$de_p<0.05,], aes_string(x="id", y="cancer", size="`-log10(P value)`", color="change")) +
    geom_point(aes())+ xlab("")+ylab("")+ 
    scale_color_manual(values = c(Hypermethylation="#CA0020",Hypomethylation="#0571B0"))+
    scale_size()+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))



a=reshape2::dcast(all_res,cancer~id,value.var = "auc")
summary(all_res$auc)
rownames(a)=a$cancer
a$cancer=NULL
a=a[,colSums(a>0.8)>1]

library(ComplexHeatmap)
library(circlize)
col_fun = circlize::colorRamp2(c(0,0.5, 1), c("#0571B0","white", "#CA0020"))
star_p=function(x){if(x>=0.8){format(x,scientific=F,digits=2)}else{" "}}
p5 = Heatmap(a, name = "AUC",column_title = "",col = col_fun,row_names_side = "left",column_names_side = "bottom", rect_gp = gpar(type = "none"),border = TRUE, 
             cell_fun = function(j, i, x, y, width, height, fill) {
                 if(a[i, j] != 0){
                     grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = col_fun(a[i, j]), fill =  col_fun(a[i, j])))
                     grid.text(star_p(a[i, j]), x, y, gp = gpar(fontsize = 8))
                 }else{grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill =  "white"))}}, cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = T, show_column_names = T)
p5

load("tcga_methy_surv_result_gene.RData")
length(unique(dat$tcga_type)) ## 33
dat=methy_surv_result[methy_surv_result$cox_p<0.05,]
dat$`-log10(COX P value)`=-log10(dat$cox_p)
dat$`-log10(COX HR)`=log10(dat$cox_hr)
p6=ggplot(dat, aes_string(y="tcga_type", x="gene", color="`-log10(COX HR)`", size="`-log10(COX P value)`")) +
    geom_point(aes())+ xlab("")+ylab("")+ scale_color_gradient(high="#CA0020",low="#0571B0")+
    scale_size()+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p6

