rm(list=ls())
options(stringsAsFactors = F)

load("medip_stad_peak_samp250_subsamp.RData")
load("traintest_norm.RData")
load("samp_select.RData")
load("medip_stad_peak_samp250_phen.RData")
phen=na.omit(phen[samp_select,c(2,3,5)])
phen$Gender=as.numeric(phen$Gender)-1

limma_de_all=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(phen)]
    coldata=phen[train_samp,]
    psm <- matchit(Group~Age+Gender,data=coldata,method="nearest",distance = "logit",ratio = 1)
    coldata <- match.data(psm)
    design=model.matrix(~0+coldata$subclass+coldata$Group)
    colnames(design)=gsub("coldata[$]","",colnames(design))
    res=lmFit(norm[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="Groupcancer") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[abs(res$logFC)>1 & res$adj.P.Val<0.01,]
    message(nrow(limma_de_all[[i]]))
}
save(limma_de_all,file="limma_PSM_de.RData")

