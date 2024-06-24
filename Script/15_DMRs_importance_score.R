load("limma_modlist_PSM_rf.RData")
all_imp=data.frame()
for (i in 1:100) {
    all_imp=rbind(all_imp,data.frame(pos=rownames(limma_modlist[[i]]$rf_mod$finalModel$importance),importance=limma_modlist[[i]]$rf_mod$finalModel$importance))
}

a=do.call(rbind,lapply(unique(all_imp$pos), function(x){data.frame(pos=x,med=median(all_imp[all_imp$pos %in% x,"MeanDecreaseGini"]))}))
a=a[order(a$med,decreasing = T),]
all_imp$pos=factor(all_imp$pos,levels = a$pos)

load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
imp_de_gene2=imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,]
imp_de_gene2$group=imp_de_gene[match(imp_de_gene2$pos,imp_de_gene$pos),"group"]

hgnc_complete_set <- read.delim("Annotation_file/hgnc_complete_set.txt")
imp_de_gene2$type=hgnc_complete_set[match(imp_de_gene2$SYMBOL,hgnc_complete_set$symbol),"locus_group"]
# CCDC144BP  Pseudogene
# LOC100134317 lncRNA
#  PTGER4P2-CDK2AP2P2 Transcribed Pseudogene
imp_de_gene2[imp_de_gene2$SYMBOL %in% c("CCDC144B","PTGER4P2-CDK2AP2P2"),"type"]="pseudogene"
imp_de_gene2[imp_de_gene2$SYMBOL %in% c("LOC100134317"),"type"]="non-coding RNA"
table(imp_de_gene2$type,imp_de_gene2$group)
all_imp=merge(all_imp,imp_de_gene2,by="pos")

p2=ggplot(all_imp, aes(pos, MeanDecreaseGini)) +
    geom_boxplot(aes(fill=group))+ geom_jitter(width = 0.1, alpha = 0.1, color = 'black')+
    scale_fill_brewer(palette = "Set1")+
    labs(title = "")+xlab(label = "")+ylab("MeanDecreaseGini")+theme(axis.text.x = element_blank())
p2
