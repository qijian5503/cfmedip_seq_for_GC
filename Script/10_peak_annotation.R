## all peak annotation
rm(list=ls())
options(stringsAsFactors = F)

load("medip_stad_peak_samp250_raw_exp.RData")
peak_ana=data.frame(pos=rownames(peak_exp),max=rowSums(peak_exp),count=rowSums(peak_exp>0))
peak_ana=peak_ana[!grepl("chrUn|random|hap",peak_ana$pos),]

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(tidyverse)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene

pos=do.call(rbind,lapply(as.character(peak_ana$pos), function(x){unlist(strsplit(x,"_"))})) %>% as.data.frame()
pos=na.omit(data.frame(seqnames=pos$V1,start=as.numeric(pos$V2),end=as.numeric(pos$V3)))
peakAnno = GRanges(seqnames=Rle(pos[,1]),ranges=IRanges(pos[,2], pos[,3]), strand=rep(c("*"), nrow(pos))) %>% annotatePeak(. ,TxDb=txdb, annoDb="org.Hs.eg.db") %>% as.data.frame()
peakAnno$annotation1=str_split(peakAnno$annotation,'[(]',simplify = T)[,1]
table(peakAnno$annotation1)

peakAnno$pos=paste(peakAnno$seqnames,peakAnno$start,peakAnno$end,sep = "_")
peak_ana=merge(peak_ana,peakAnno[,c(6,14:16,18:19)],by="pos")
save(peak_ana,file="medip_stad_peakAnno.RData")

table(peak_ana$annotation1)
a=as.data.frame(table(peak_ana$annotation1))
colnames(a)=c("Region","Freq")
a$Freq=a$Freq/nrow(peak_ana)*100

ggplot(data = a,aes(x=Region,y=Freq,fill=Region))+
    geom_bar(stat="identity")+xlab(label = "")+ylab("Abundance")



## genomic repeats regions annotation
library(annotatr)
repeat_ucsc_table_hg19 <- read.delim("Annotation_file/repeat_ucsc_table_hg19.csv", header=FALSE, comment.char="#")
repeat_ucsc_table_hg19=repeat_ucsc_table_hg19[,c(6:8,11:13)]
colnames(repeat_ucsc_table_hg19)=NULL
write.table(repeat_ucsc_table_hg19,"repeat_ucsc_table_hg19.bed",quote = F,sep = "\t",row.names = F,col.names = F)

load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
imp_de_gene2=imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,]
imp_de_gene=imp_de_gene[imp_de_gene$pos %in% imp_de_gene2$pos,]
imp_de_gene=unique(imp_de_gene[,c(1,5,6,7)])
de_gene=as.data.frame(str_split(imp_de_gene$pos,"_",simplify = T))
gpl = GRanges(seqnames=Rle(de_gene[,1]),ranges=IRanges(as.numeric(de_gene[,2]), as.numeric(de_gene[,3])), strand=rep(c("*"), nrow(de_gene))) 

extraCols = c(repName = 'character', repClass = 'character', repFamily = 'character')
repeat_regions = read_regions(con = "repeat_ucsc_table_hg19.bed", genome = 'hg19', extraCols = extraCols, format = 'bed')
repeat_annotated = annotate_regions(regions = gpl,annotations = repeat_regions,ignore.strand = TRUE,quiet = FALSE)
repeat_annotated=repeat_annotated %>% as.data.frame()
repeat_annotated=unique(repeat_annotated[,c(1:3,11:13)])
repeat_annotated$pos=paste(repeat_annotated$seqnames,repeat_annotated$start,repeat_annotated$end,sep = "_")
table(repeat_annotated$pos %in% imp_de_gene$pos)
imp_de_gene=merge(imp_de_gene,repeat_annotated,by="pos",all=T)
table(unique(imp_de_gene[imp_de_gene$group %in% "Hypomethylation","group"]))
save(imp_de_gene,file="medip_stad_repeat_annotated.RData")


table(imp_de_gene$group,imp_de_gene$annot.repClass)
p_value=chisq.test(table(imp_de_gene$group,imp_de_gene$annot.repClass))$p.value
a= imp_de_gene[!is.na(imp_de_gene$annot.repClass),c("annot.repClass","group")] %>% dplyr::group_by(group,annot.repClass) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
a$n_lab=ifelse(a$Freq<1,"",a$n)
a$Freq_lab=ifelse(a$Freq<1,"",paste(round(a$Freq,2),"%"))
ggplot(a,aes(x=group,y=Freq,fill=annot.repClass))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))+
    labs(subtitle = paste0("p = ",round(p_value,3)),title="", x="", y="Frenquency of response") 


## 58 DMRs annotation
load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
imp_de_gene2=imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,]
imp_de_gene2$group=imp_de_gene[match(imp_de_gene2$pos,imp_de_gene$pos),"group"]

hgnc_complete_set <- read.delim("hgnc_complete_set.txt")
table(imp_de_gene2$SYMBOL %in% hgnc_complete_set$symbol)
imp_de_gene2$SYMBOL[!imp_de_gene2$SYMBOL %in% hgnc_complete_set$symbol ]

imp_de_gene2$type=hgnc_complete_set[match(imp_de_gene2$SYMBOL,hgnc_complete_set$symbol),"locus_group"]
imp_de_gene2[imp_de_gene2$SYMBOL %in% c("CCDC144B","LOC653513","LOC441666","PTGER4P2-CDK2AP2P2"),"type"]="pseudogene"
imp_de_gene2[imp_de_gene2$SYMBOL %in% c("LOC100134317"),"type"]="non-coding RNA"

a= imp_de_gene2[,c("type","group")] %>% dplyr::group_by(group,type) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
a$n_lab=ifelse(a$Freq<1,"",a$n)
a$Freq_lab=ifelse(a$Freq<1,"",paste(round(a$Freq,2),"%"))
a$group=gsub("methylation","",a$group)
p2=ggplot(a,aes(x=group,y=Freq,fill=type))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))


b= imp_de_gene2[,c("group","annotation1")] %>% dplyr::group_by(group,annotation1) %>% 
    dplyr::count() %>% group_by(group) %>% mutate(Freq=n/sum(n)*100)
b$n_lab=ifelse(b$Freq<1,"",b$n)
b$Freq_lab=ifelse(b$Freq<1,"",paste(round(b$Freq,2),"%"))
b$group=gsub("methylation","",b$group)

p3=ggplot(b,aes(x=group,y=Freq,fill=annotation1))+geom_col()+geom_text(aes(label=Freq_lab),position=position_stack(vjust=0.5))+theme_classic()+
    scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"))+
    labs(x = "",y = "Frequency")+
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0))
