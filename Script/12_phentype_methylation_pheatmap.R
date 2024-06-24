library(pheatmap)
library(gplots)
library(colorRamps)
library(RColorBrewer)

load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
table(imp_de_gene2$Hypermethylation>10)
table(imp_de_gene2$Hypomethylation>10)
de_pos=as.character(imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,"pos"])

load("traintest_norm.RData")
load("valid_norm.RData")
norm=cbind(norm,norm_v)
nmf_exp=norm[de_pos,]

load("medip_stad_peak_samp250_phen.RData")
table(phen$sample %in% colnames(nmf_exp))
phen=phen[phen$sample %in% colnames(nmf_exp),]
colnames(phen)
phen=phen[,c(5,6,3,2,8,9)]
phen=phen[order(phen$Group),]

phen$Group=factor(phen$Group,labels = c("Normal","Cancer"),levels = c("normal","cancer"))
phen$Stage=factor(phen$Stage,labels = c("Normal","Stage 0","Stage I","Stage II","Stage III","Stage IV"),levels = c("normal","0","I","II","III","IV"))
median(phen$Age,na.rm = T)
phen$Age=ifelse(phen$Age>=56,">=56",ifelse(phen$Age<56,"<56",NA))
phen$Age=factor(phen$Age,levels = c("<56",">=56"))
phen$Subtype=factor(phen$Subtype)
phen$Location=factor(phen$Location)


ann_colors = list(
    Group=c(Cancer="#CA0020", Normal="#0571B0"),
    change=c(Hypermethylation="#CA0020", Hypomethylation="#0571B0"),
    Gender=c(`F`="#FB9A99",`M`="#A65628"),
    Age=c(`>=56`="#4DAF4A",`<56`="#984EA3"),
    Stage=c(`Normal`="#0571B0",`Stage 0`="#F5F5F5",`Stage I`="#C7EAE5",`Stage II`="#80CDC1",`Stage III`="#35978F",`Stage IV`="#01665E",`Stage V`="#003C30"),
    Subtype=c(DGC="#00468BFF", IGC="#42B540FF", Mix="#0099B4FF"),
    Location=c(antrum="#CA0020",body="#00468BFF",cardia="#42B540FF",`cardia to body`="#0099B4FF",`cardia to fundus`="#FDAF91FF")
)
phen1=data.frame(row.names = rownames(phen),Group=phen$Group)

split_row=unique(imp_de_gene[c(1,5)])
rownames(split_row)=split_row$pos
split_row=split_row[de_pos,]
colnames(split_row)[2]="change"
split_row$pos=NULL

p5=ComplexHeatmap::pheatmap(as.matrix(nmf_exp[,rownames(phen)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors,
                            col = colorRampPalette(c("navy","white","firebrick3"))(10),
                            annotation_col = phen,annotation_row = split_row,
                            row_split = split_row$change,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)
p5


