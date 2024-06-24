## AUC Cohort II       ##########################################
rm(list=ls())

load("limma_modlist_PSM_rf.RData")
load("medip_stad_peak_samp250_subsamp.RData")

all_test_pre=data.frame()
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    all_test_pre=rbind(all_test_pre,limma_pre[limma_pre$sample %in% test_samp,])
    all_valid_pre=rbind(all_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}


library(rms)
rcorrcens(all_valid_pre$group~all_valid_pre$pre)
#                   C    Dxy  aDxy    SD      Z P    n
# all_valid_pre$pre 0.98 0.96 0.96 0.003 322.41 0 4600

ci(roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre)))
# 95% CI: 0.9773-0.9832 (DeLong)

library(pROC)
coords(roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre)), "best", ret = "all", transpose = FALSE)
# 0.527     0.9423333   0.88375      0.9219565

g2=pROC::ggroc(list(`Validation set`=roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre))))
auc2=roc(as.factor(all_valid_pre$group), as.numeric(all_valid_pre$pre))


p1=g2+ggsci::scale_color_lancet()+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    theme_classic() + theme(legend.position = "")+
    annotate("text", y = 0.5, x=0.5, label = paste0("Cohort II\n AUC = ",round(auc2$auc[[1]],4),"\n95% CI ( ",round(ci(auc2)[1],4)," - ",round(ci(auc2)[3],4)," )"),colour ="#00468BFF")



## phen risk roc valid ########################
load("limma_modlist_PSM_rf.RData")
load("medip_stad_peak_samp250_subsamp.RData")
all_test_pre=data.frame()
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    all_test_pre=rbind(all_test_pre,limma_pre[limma_pre$sample %in% test_samp,])
    all_valid_pre=rbind(all_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}
all_pre=do.call(rbind,lapply(unique(all_valid_pre$sample),function(x){
    data.frame(sample=x,risk=median(all_valid_pre[all_valid_pre$sample %in% x,"pre"]))
}))

load("medip_stad_peak_samp250_phen.RData")
table(phen$Stage)
phen=merge(phen,all_pre,by="sample")

for (i in c("III","I","IV","II","0")) {
    dat=phen[phen$Stage %in% c("normal",i),]
    assign(paste0("Stage",i,"_auc"),roc(as.factor(dat$Group), as.numeric(dat$risk)))
    assign(paste0("Stage",i,"_auc1"),roc(as.factor(dat$Group), as.numeric(dat$risk))$auc[[1]])
}
g1=pROC::ggroc(list(`Stage 0`=Stage0_auc,`Stage I`=StageI_auc,`Stage II`=StageII_auc,
                    `Stage III`=StageIII_auc,`Stage IV`=StageIV_auc))
p2=g1+ggsci::scale_color_lancet() + labs(x = "1 - Specificity", y = "Sensitivity") + 
    theme_classic() + theme(legend.title = element_blank())+
    annotate("text", y = 0.5, x=0.5, 
             label = paste0("Cohort II\nStage 0 AUC = ",round(Stage0_auc1,4)," \nStage I AUC = ",round(StageI_auc1,4)," \nStage II AUC = ",round(StageII_auc1,4)," \nStage III AUC = ",round(StageIII_auc1,4)," \nStage IV AUC = ",round(StageIV_auc1,4)),colour ="#00468B99")
p2



g1=pROC::ggroc(list(Age=roc(as.factor(phen$Group), as.numeric(phen$Age)),
                    Gender=roc(as.factor(phen$Group), as.numeric(phen$Gender)),
                    Concentration=roc(as.factor(phen$Group), as.numeric(phen$Concentration)),
                    Risk=roc(as.factor(phen$Group), as.numeric(phen$risk))))
auc1=roc(as.factor(phen$Group), as.numeric(phen$Age))$auc[[1]]
auc2=roc(as.factor(phen$Group), as.numeric(phen$Gender))$auc[[1]]
auc3=roc(as.factor(phen$Group), as.numeric(phen$Concentration))$auc[[1]]
auc4=roc(as.factor(phen$Group), as.numeric(phen$risk))$auc[[1]]

p3=g1+ggsci::scale_color_lancet() + labs(x = "1 - Specificity", y = "Sensitivity") + 
    theme_classic() + theme(legend.title = element_blank())+
    annotate("text", y = 0.5, x=0.5, 
             label = paste0("Cohort II\n                 Age AUC = ",round(auc1,4)," \n           Gender AUC = ",round(auc2,4)," \nConcentration AUC = ",round(auc3,4)," \n               Risk AUC = ",round(auc4,4)),colour ="#00468B99")
p3



## sub-group valid     ######################################
load("limma_modlist_PSM_rf.RData")
load("medip_stad_peak_samp250_subsamp.RData")
all_valid_pre=data.frame()
for (i in 1:100) {
    limma_pre=limma_modlist[[i]]$limma_pre
    test_samp=subsample[[i]]$test_samp
    all_valid_pre=rbind(all_valid_pre,limma_pre[!limma_pre$sample %in% test_samp,])
}
all_pre=do.call(rbind,lapply(unique(all_valid_pre$sample),function(x){
    data.frame(sample=x,risk=median(all_valid_pre[all_valid_pre$sample %in% x,"pre"]))
}))

load("medip_stad_peak_samp250_phen.RData")
table(phen$Stage)
phen=merge(phen,all_pre,by="sample")

median(phen$Age,na.rm = T)
phen$Age=ifelse(phen$Age<56,"<56",ifelse(phen$Age>=56,">=56",NA))
median(phen$Concentration,na.rm = T)
phen$Concentration=ifelse(phen$Concentration<7.25,"<7.25",ifelse(phen$Concentration>=7.25,">=7.25",NA))
#phen[is.na(phen$Concentration),"Concentration"]="Unknown"
#phen[is.na(phen$Age),"Age"]="Unknown"
#phen$Gender=as.character(phen$Gender)
#phen[is.na(phen$Gender),"Gender"]="Unknown"

dat=reshape2::melt(phen[,c(1:5,8)],c("sample","risk","Group"))
dat$sub_group=paste0(dat$variable," ",dat$value)
dat$risk=1-dat$risk
dat=dat[!is.na(dat$value),]

dat2=reshape2::melt(phen[,c(1,5,6,8)],c("sample","risk","Group"))
dat2=dat2[!is.na(dat2$value),]
dat2$sub_group=ifelse(!dat2$value %in% "normal",paste0(dat2$variable," ",dat2$value),dat2$value)
dat2$risk=1-dat2$risk

dat3=reshape2::melt(phen[,c(1,5,6,8)],c("sample","risk","Group"))
dat3=dat3[!is.na(dat3$value),]
dat3=dat3[!dat3$value %in% "normal",]
dat3$sub_group=ifelse(!dat3$value %in% "normal",paste0(dat3$variable," ",dat3$value),dat3$value)
dat3$risk=1-dat3$risk


p4=ggplot(dat, aes(sub_group, risk,fill=sub_group)) + 
    geom_boxplot()+geom_jitter(width = 0.1, alpha = 0.5, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    stat_compare_means(aes(group=sub_group),label.y = 1.1,label = "p.format")+
    facet_grid(.~Group)+theme_bw()+ylab("Predicted risk score")+
    #theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),panel.grid = element_blank())
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),legend.title = element_blank())
p4
p5=ggplot(dat3, aes(sub_group, risk,fill=sub_group)) + 
    geom_boxplot()+geom_jitter(width = 0.1, alpha = 0.5, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    stat_compare_means(aes(group=sub_group),label.y = 1.1,label = "p.format")+
    theme_bw()+ylab("Predicted risk score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),legend.title = element_blank())
p5


