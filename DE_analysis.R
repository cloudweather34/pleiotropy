setwd("/Volumes/cluster/Wei-Yun/project2/WGCNA/")
network=read.delim("/Volumes/cluster/Wei-Yun/network/fly_genetic_interactions.txt",sep = "\t")
rm(list=ls())
library(edgeR)
library(pheatmap)
library(biomaRt)
library(pathview)
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(gplots)
library(magrittr)
library(dplyr)
library(ExpressionNormalizationWorkflow)
library(car)
library(colorspace)
library(readxl)
library(scales)
library(gtools)
library(pheatmap)
library(topGO)
library(GOSim)

ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}

####replicate-specific contrasts#### (DE genes)####
dat_B_H=read.csv(file = "~/Dropbox (PopGen)/backup/Wei-Yun/monster_fullset_readcounts.csv",sep = ";")
colnames(dat_B_H)[1]="name"
dat_B_H_filtered=dat_B_H[apply(cpm(dat_B_H[,-1]),1,function(x){!any(x<=1)}),]
group=c("H1.1","H2.1","H3.1","H4.1","H5.1","H6.1","H7.1","H8.1","H9.1","H10.1","B27","B30",
        "H1.2","H2.2","H3.2","H4.2","H5.2","H6.2","H7.2","H8.2","H9.2","H10.2","B28","H1.3",
        "H2.3","H3.3","H4.3","H5.3","H6.3","H7.3","H9.3","H10.3","B26","B29","H8.3")

row.names(dat_B_H_filtered)=dat_B_H_filtered[,1]
dat_B_H_filtered=dat_B_H_filtered[,-1]
dat_B_H_filtered_log=log(dat_B_H_filtered)
group_rep=paste(substr(colnames(dat_B_H_filtered),1,1),substr(colnames(dat_B_H_filtered),8,8),sep = "_")
group_rep[c(11,12,23,33,34)]="B"
y=DGEList(counts = dat_B_H_filtered,group = group_rep)
ModelDesign_rep=model.matrix(~0+group_rep)
DGE_rep=estimateDisp(y,design = ModelDesign_rep,robust = T)
GLM_rep=glmFit(DGE_rep,design = ModelDesign_rep)
res_table=list()
sig_ID_all=list()
sig_ID_up=list()
sig_ID_dn=list()

for (i in 1:10){
  my.contrasts_m=rep(0,length(unique(group_rep)))
  my.contrasts_m[c(1,i+1)]=c(-1,1)
  LRT_m=glmLRT(GLM_rep,contrast = my.contrasts_m)
  res_m=LRT_m$table
  row.names(res_m)=row.names(dat_B_H_filtered)
  res_m$padj=p.adjust(res_m$PValue,method = "BH")
  res_table[[i]]=res_m
  sig_ID_up[[i]]=row.names(res_m[res_m$padj<0.05&res_m$logFC>0,])
  sig_ID_dn[[i]]=row.names(res_m[res_m$padj<0.05&res_m$logFC<0,])
  sig_ID_all[[i]]=row.names(res_m[res_m$padj<0.05,])
#  write.table(res_m[res_m$padj<0.05,],file = paste0("/Users/weiyun/Dropbox (PopGen)/Wei-Yun (1)/manuscript_targetofselection/table/DE results/DE_",i,".txt"),quote = F,sep = "\t")
}

setwd("/Volumes/cluster/Wei-Yun/metabolomic_data/")


####jaccard index of DE genes####
ind=combinations(10,2,set=TRUE, repeats.allowed=FALSE)
ja_de=c()
for (i in 1:45) {
  inte=length(intersect(unlist(sig_ID_all[ind[i,1]]),unlist(sig_ID_all[ind[i,2]])))
  uni=length(union(unlist(sig_ID_all[ind[i,1]]),unlist(sig_ID_all[ind[i,2]])))
  ja_de[i]=inte/uni
}


####replicate-specific GO analysis####
setwd("/Volumes/cluster/Wei-Yun/project2/GOdatabase/")
GO_res_table_up=list()
for (i in 1:10){
  tmp=factor(as.integer(rownames(dat_B_H_filtered)%in%unlist(sig_ID_up[[i]])))
  names(tmp)=rownames(dat_B_H_filtered)
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.classic",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  #  tmp_res$Fisher.classic=resTopGO.classic@score
  GO_res_table_up[[i]]=tmp_res
}


GO_res_table_up=lapply(GO_res_table_up,function(x) {
  x$Fisher.classic[x$Fisher.classic=="< 1e-30"]=1e-30
  x$Fisher.classic.padj=p.adjust(x$Fisher.classic,method = "BH")
  return(x)})

GO_res_table_dn=list()
for (i in 1:10){
  tmp=factor(as.integer(rownames(dat_B_H_filtered)%in%unlist(sig_ID_dn[[i]])))
  names(tmp)=rownames(dat_B_H_filtered)
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.classic",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  #  tmp_res$Fisher.classic=resTopGO.classic@score
  GO_res_table_dn[[i]]=tmp_res
}

GO_res_table_dn=lapply(GO_res_table_dn,function(x) {
  x$Fisher.classic[x$Fisher.classic=="< 1e-30"]=1e-30
  x$Fisher.classic.padj=p.adjust(x$Fisher.classic,method = "BH")
  return(x)})

####Jaccard Index of enriched up/down-regulated GO terms####
ind=combinations(10,2,set=TRUE, repeats.allowed=FALSE)
GO.use.down=lapply(GO_res_table_dn,function(x) x$GO.ID[which(x$Fisher.classic.padj<0.05)])
GO.use.up=lapply(GO_res_table_up,function(x) x$GO.ID[which(x$Fisher.classic.padj<0.05)])

ja_goup=c()
for (i in 1:45) {
  inte=length(intersect(unlist(GO.use.up[ind[i,1]]),unlist(GO.use.up[ind[i,2]])))
  uni=length(union(unlist(GO.use.up[ind[i,1]]),unlist(GO.use.up[ind[i,2]])))
  ja_goup[i]=inte/uni
}

ja_godn=c()
for (i in 1:45) {
  inte=length(intersect(unlist(GO.use.down[ind[i,1]]),unlist(GO.use.down[ind[i,2]])))
  uni=length(union(unlist(GO.use.down[ind[i,1]]),unlist(GO.use.down[ind[i,2]])))
  ja_godn[i]=inte/uni
}

####metabolomic analysis####
setwd("/Volumes/cluster/Wei-Yun/metabolomic_data/")
dat.met=read.csv(file = "metabolome.csv")
dat.met.bhqc=dat.met[,c(1:8,24:41,50:57)]
colnames(dat.met.bhqc)[-c(1:3)]=c(paste0("B",1:5),paste0("H",rep(4:9,each=3),rep(1:3,6)),paste0("QC",1:8))

met.test=list()
ind=seq(6,23,3)
dat.met.bhqc.use=log(dat.met.bhqc[,-c(1:3)])

for (i in 1:6) {
  temp=matrix(NA,nrow = 940,ncol = 4)
  for (j in 1:940) {
    a=t.test(dat.met.bhqc.use[j,1:5],dat.met.bhqc.use[j,c(ind[i]:(ind[i]+2))])
    temp[j,1]=a$estimate[2]
    temp[j,2]=a$statistic
    temp[j,3]=a$p.value
  }
  temp[,4]=p.adjust(temp[,3],method = "BH")
  met.test[[i]]=temp
}


ja.use=lapply(met.test, function(x) dat.met.bhqc[x[,4]<0.1,1])
sapply(ja.use,length)


library(gtools)
ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE)
ja_met=c()
for (j in 1:15) {
  inte=length(intersect(unlist(ja.use[ind[j,1]]),unlist(ja.use[ind[j,2]])))
  uni=length(union(unlist(ja.use[ind[j,1]]),unlist(ja.use[ind[j,2]])))
  ja_met[j]=inte/uni
}

####gene expression data - hot-evolved####
dat_ge_hotevo=dat_B_H_filtered[,c(1:10,13:22,24:32,35)]
dat_cor1=cor(dat_ge_hotevo[,c(1,11,21)])[upper.tri(cor(dat_ge_hotevo[,c(1,11,21)]))]
dat_cor2=cor(dat_ge_hotevo[,c(2,12,22)])[upper.tri(cor(dat_ge_hotevo[,c(2,12,22)]))]
dat_cor3=cor(dat_ge_hotevo[,c(3,13,23)])[upper.tri(cor(dat_ge_hotevo[,c(3,13,23)]))]
dat_cor4=cor(dat_ge_hotevo[,c(4,14,24)])[upper.tri(cor(dat_ge_hotevo[,c(4,14,24)]))]
dat_cor5=cor(dat_ge_hotevo[,c(5,15,25)])[upper.tri(cor(dat_ge_hotevo[,c(5,15,25)]))]
dat_cor6=cor(dat_ge_hotevo[,c(6,16,26)])[upper.tri(cor(dat_ge_hotevo[,c(6,16,26)]))]
dat_cor7=cor(dat_ge_hotevo[,c(7,17,27)])[upper.tri(cor(dat_ge_hotevo[,c(7,17,27)]))]
dat_cor8=cor(dat_ge_hotevo[,c(8,18,30)])[upper.tri(cor(dat_ge_hotevo[,c(8,18,30)]))]
dat_cor9=cor(dat_ge_hotevo[,c(9,19,28)])[upper.tri(cor(dat_ge_hotevo[,c(4,14,28)]))]
dat_cor10=cor(dat_ge_hotevo[,c(10,20,29)])[upper.tri(cor(dat_ge_hotevo[,c(10,20,29)]))]
cor_gen_within=c(dat_cor1,dat_cor2,dat_cor3,dat_cor4,dat_cor5,dat_cor6,dat_cor7,dat_cor8,dat_cor9,dat_cor10)

dat.evo.mean=data.frame(evo.1=apply(dat_ge_hotevo[,c(1,11,21)],1,mean),
                        evo.2=apply(dat_ge_hotevo[,c(2,12,22)],1,mean),
                        evo.3=apply(dat_ge_hotevo[,c(3,13,23)],1,mean),
                        evo.4=apply(dat_ge_hotevo[,c(4,14,24)],1,mean),
                        evo.5=apply(dat_ge_hotevo[,c(5,15,25)],1,mean),
                        evo.6=apply(dat_ge_hotevo[,c(6,16,26)],1,mean),
                        evo.7=apply(dat_ge_hotevo[,c(7,17,27)],1,mean),
                        evo.8=apply(dat_ge_hotevo[,c(8,18,30)],1,mean),
                        evo.9=apply(dat_ge_hotevo[,c(9,19,28)],1,mean),
                        evo.10=apply(dat_ge_hotevo[,c(10,20,29)],1,mean))

dat_gen_cor_across=cor(dat.evo.mean)[upper.tri(cor(dat.evo.mean))]

boxplot(dat_gen_cor_across,cor_gen_within)

library(biomaRt)
ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}

conv_query=lapply(sig_ID_all,function(x) ID_converter(ID = x,db = ensembl,attributes = c("external_gene_name","entrezgene_id"),filters = "flybase_gene_id")[,2])

all_bg=row.names(res_table[[1]])
conv_query_bg=ID_converter(ID = all_bg,db = ensembl,attributes = c("flybase_gene_id","entrezgene_id"),filters = "flybase_gene_id")


for (i in 1:10) {
  write.table(conv_query[[i]],file = paste0(i,"_gene.candidate.entrezgene.txt"),
              quote = F,row.names = F)
}

res_table_output=list()

for (i in 1:10) {
  temp=data.frame(row.names(res_table[[i]]),res_table[[i]][,1])
  colnames(temp)=c("flybase_gene_id","logFC")
  res_table_output[[i]]=merge(conv_query_bg,temp)[,c(2,3)]
  write.table(res_table_output[[i]],file = paste0(i,"_gene.logfc.candidate.txt"),
              quote = F,row.names = F,col.names = F)
  }


engry_meta=genesInTerm(tgd1,"GO:0009150")[[1]]
ion_transport=genesInTerm(tgd1,"GO:0098662")[[1]]

group_rep_evo=group_rep[-which(group_rep=="B")]
ind=unique(group_rep_evo)
dat.avg.fc=matrix(NA,10780,10)

for (i in 1:10) {
  dat.avg.fc[,i]=apply(dat_fc[,which(group_rep_evo==ind[i])],1,mean)
}

row.names(dat.avg.fc)=row.names(dat_fc)

bg_cor=cor(dat.avg.fc)[upper.tri(cor(dat.avg.fc))]
en_cor=cor(dat.avg.fc[engry_meta,])[upper.tri(cor(dat.avg.fc[engry_meta,]))]
ion_cor=cor(dat.avg.fc[ion_transport,])[upper.tri(cor(dat.avg.fc[ion_transport,]))]

boxplot(en_cor,ion_cor)



png("/Volumes/cluster/Wei-Yun/project2/hetero_GO.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
plot_dat=data.frame(h=c(1-bg_cor,
                        1-ion_cor,
                        1-en_cor),
                    cat=as.factor(c(rep("Background",45),rep("Ion transport",45),rep("Purine metabolism",45))))
qplot(cat,h,data = plot_dat,geom = "violin",fill=cat,alpha=I(0.5),xlab = "",ylab="Heterogeneity across replicates")+
  theme_classic()+
  theme(axis.text.x=element_text(size=8),legend.position = "NULL")+
  geom_boxplot(width=0.1,fill="white")+
  scale_fill_manual(values=c("grey","salmon","cyan"))
dev.off()