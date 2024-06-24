
load("imp_de_gene.RData")
imp_de_gene2=reshape2::dcast(imp_de_gene,pos+SYMBOL+annotation1~group)
imp_de_gene2=imp_de_gene2[imp_de_gene2$Hypermethylation>10 | imp_de_gene2$Hypomethylation>10,]
imp_de_gene=imp_de_gene[imp_de_gene$pos %in% imp_de_gene2$pos,]
imp_de_gene=unique(imp_de_gene[,c(1,5,6,7)])
table(imp_de_gene$group)

down_gene=unique(imp_de_gene[imp_de_gene$group %in% "Hypomethylation","SYMBOL"])
up_gene=unique(imp_de_gene[imp_de_gene$group %in% "Hypermethylation","SYMBOL"])

library(clusterProfiler)
library(org.Hs.eg.db)
ensem_down=bitr(down_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
go_down=enrichGO(ensem_down$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "none",pvalueCutoff = 0.1,qvalueCutoff = 0.2,keyType = "ENTREZID")
kegg_down=enrichKEGG(ensem_down$ENTREZID, organism = 'hsa', keyType = 'kegg',pAdjustMethod = "none", pvalueCutoff = 0.1,qvalueCutoff = 0.2)
dotplot(kegg_down, showCategory=15) 
kegg_down = setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID") ## 显示基因symbol

ensem_up=bitr(up_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
go_up=enrichGO(ensem_up$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "none",pvalueCutoff = 0.1,qvalueCutoff = 0.2,keyType = "ENTREZID")
kegg_up=enrichKEGG(ensem_up$ENTREZID, organism = 'hsa', keyType = 'kegg',pAdjustMethod = "none", pvalueCutoff = 0.1,qvalueCutoff = 0.2)
dotplot(kegg_up, showCategory=15) 
kegg_up = setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID") ## 显示基因symbol

enrich_res=list(go_down=go_down,kegg_down=kegg_down,go_up=go_up,kegg_up=kegg_up)
save(enrich_res,file = "enrich_up_down_res58.RData")


load("enrich_up_down_res58.RData")
go_up=enrich_res$go_up
kegg_up=enrich_res$kegg_up
go_down=enrich_res$go_down
kegg_down=enrich_res$kegg_down

p1=dotplot(kegg_up, showCategory=10) 
p2=dotplot(go_up, split="ONTOLOGY", showCategory=3)+facet_grid(ONTOLOGY~., scale="free")

#dotplot(kegg_down, showCategory=10) 
p3=dotplot(go_down, split="ONTOLOGY", showCategory=3)+facet_grid(ONTOLOGY~., scale="free")
