## subtype             ########################################
load("samp_select.RData")
load("medip_stad_peak_samp250_phen.RData")
unique(phen$Location)
my_comparisons=list(c("antrum","body"),c("antrum","cardia"),c("cardia","body"),
                    c("cardia to fundus","cardia to body"),c("cardia","cardia to body"))
ggplot(phen[!phen$Location %in% c(NA,""),], aes(Location, Concentration)) + 
    geom_boxplot(aes(fill=Location))+xlab(label = "")+
    stat_compare_means(comparisons = my_comparisons,aes(group=Location),hide.ns = T)+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))

unique(phen$Subtype)
my_comparisons2=list(c("DGC","IGC"),c("DGC","Mix"),c("IGC","Mix"))
ggplot(phen[!phen$Subtype %in% c(NA,""),], aes(Subtype, Concentration)) + 
    geom_boxplot(aes(fill=Subtype))+xlab(label = "")+
    scale_fill_manual(values=c("#0571B0","#CA0020","#42B540FF"))+
    stat_compare_means(comparisons = my_comparisons2,aes(group=Subtype),hide.ns = T)+labs(title = "")+ylab("cfDNA Concentration(ng/ml)")+
    geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))



load("traintest_norm.RData")
load("valid_norm.RData")
norm=norm[,colnames(norm) %in% phen$sample]
phen=phen[phen$sample %in% colnames(norm),]

coldata=phen[!phen$Subtype %in% c(NA,""),]
design=model.matrix(~0+coldata$Subtype)
colnames(design)=gsub("coldata[$]Subtype","",colnames(design))
res1=lmFit(norm[,coldata$sample],design = design,) %>% contrasts.fit(., makeContrasts(DGC-IGC, levels = design)) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
res2=lmFit(norm[,coldata$sample],design = design,) %>% contrasts.fit(., makeContrasts(DGC-Mix, levels = design)) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
res3=lmFit(norm[,coldata$sample],design = design,) %>% contrasts.fit(., makeContrasts(IGC-Mix, levels = design)) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
table(abs(res1$logFC)>1 & res1$adj.P.Val<0.001)
table(abs(res2$logFC)>1 & res2$adj.P.Val<0.001)
table(abs(res3$logFC)>1 & res3$adj.P.Val<0.001)
# FALSE  TRUE 
# 30730     3 
save(res1,file = "subtype_de_DGC_IGC.RData")
save(res2,file = "subtype_de_DGC_Mix.RData")
save(res3,file = "subtype_de_IGC_Mix.RData")


