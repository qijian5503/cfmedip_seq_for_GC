load("Data/limma_modlist_PSM_rf.RData")
load("Data/medip_stad_peakAnno.RData")
load("Data/limma_PSM_de.RData")
limma_de=list()
for(i in 1:100){
    res=limma_de_all[[i]][1:1000,]
    limma_de=rbind(limma_de,data.frame(res,pos=rownames(res)))
}
rownames(limma_de)=NULL

limma_de$group=ifelse(limma_de$logFC>0,"Hypermethylation","Hypomethylation")
limma_de=merge(limma_de,peak_ana[,c(1,7,8)],by="pos",all.x=T)
save(limma_de,file = "Data/limma_de_top1000.RData")

limma_de2=reshape2::dcast(limma_de,pos+SYMBOL+annotation1~group)
table(limma_de2$Hypomethylation>0)
table(limma_de2$Hypomethylation>0,limma_de2$Hypermethylation>0)
table(limma_de2$Hypomethylation>10,limma_de2$Hypermethylation>10)
limma_de2$group=limma_de[match(limma_de2$pos,limma_de$pos),"group"]
limma_de2=limma_de2[limma_de2$Hypermethylation>10 | limma_de2$Hypomethylation>10,]

hgnc_complete_set <- read.delim("Annotation_file/hgnc_complete_set.txt")
table(limma_de2$SYMBOL %in% hgnc_complete_set$symbol)
g=unique(limma_de2$SYMBOL[!limma_de2$SYMBOL %in% hgnc_complete_set$symbol ])

load("Annotation_file/gencode_annotation_gtf22.RData")
unique(g[!g %in% gtf22$gene_name ])

limma_de2$type=hgnc_complete_set[match(limma_de2$SYMBOL,hgnc_complete_set$symbol),"locus_group"]
limma_de2$type=ifelse( !is.na(limma_de2$type),limma_de2$type,gtf22[match(limma_de2$SYMBOL,gtf22$gene_name),"gene_type"])
g=unique(limma_de2[is.na(limma_de2$type),"SYMBOL"])

limma_de2[limma_de2$SYMBOL %in% c("CCDC144B","LOC100133920" ,"LOC286297","LOC650226" ,"LOC407835","DKFZP586I1420","LOC644669","LOC390705","LOC613038","LOC100129138","LOC653513","LOC441666","PTGER4P2-CDK2AP2P2"),"type"]="pseudogene"
limma_de2[limma_de2$SYMBOL %in% c("LOC100134317","LOC643441"),"type"]="non-coding RNA"
limma_de2[limma_de2$SYMBOL %in% c("LOC643441","LOC100130298","LOC389641","LOC339975","LOC284412","TSG1","LOC339298","LOC100507351"),"type"]="lncRNA"
limma_de2[is.na(limma_de2$SYMBOL) ,"type"]="other"
unique(limma_de2$type)
limma_de2$type=ifelse(limma_de2$type %in% c("protein-coding gene","protein_coding"),"protein-coding gene",
                      ifelse(limma_de2$type %in% c("lncRNA","lincRNA","non-coding RNA"),"non-coding RNA",
                             ifelse(limma_de2$type %in% c("pseudogene","unprocessed_pseudogene"),"pseudogene",limma_de2$type)))

library(tidyverse)

a= limma_de2[,c("type","group")] %>% dplyr::group_by(group,type) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
b= limma_de2[,c("group","annotation1")] %>% dplyr::group_by(group,annotation1) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
a$n_lab=ifelse(a$Freq<5,"",a$n)
a$Freq_lab=ifelse(a$Freq<5,"",paste(round(a$Freq,2),"%"))
b$n_lab=ifelse(b$Freq<5,"",b$n)
b$Freq_lab=ifelse(b$Freq<5,"",paste(round(b$Freq,2),"%"))
a$group=gsub("methylation","",a$group)
b$group=gsub("methylation","",b$group)

p9=ggplot(a[!is.na(a$type),],aes(x=group,y=Freq,fill=type))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))
p10=ggplot(b,aes(x=group,y=Freq,fill=annotation1))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))


