setwd()
BiocManager::install("gage")
library(gage)

x=strsplit(readLines("/Volumes/cluster/Wei-Yun/database/dme00001.keggList_fbgnID.txt"), "[[:space:]]+")
kegg.gs=lapply(x, tail, n=-1) 
names(kegg.gs)=lapply(x, head, n=1)
rm(x)

ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
res.ID=getBM(attributes=c('flybase_gene_id','ensembl_gene_id','entrezgene_id'),filters='flybase_gene_id',mart=ensembl,values=row.names(dat_B_H_filtered))

setwd("/Volumes/cluster/Wei-Yun/project2/gage_output/")
require(gage)

gage_up=list()
gage_up_sig=list()
for (i in 1:10) {
  gage_up[[i]]=gage(cpm(y),gsets=kegg.gs,ref=which(group=="B"),samp=which(group==paste0("H",i)),rank.test = T,compare = "unpaired")
  gage_up_sig[[i]]=sigGeneSet(gage_up[[i]],outname=paste0("H",i,".kegg"),cutoff=0.05,heatmap = T)
}

gage_dn=list()
gage_dn_sig=list()
for (i in 1:10) {
  gage_dn[[i]]=gage(cpm(y),gsets=kegg.gs,ref=which(group=="B"),samp=which(group==paste0("H",i)),rank.test = T,compare = "unpaired")
  gage_dn_sig[[i]]=sigGeneSet(gage_up[[i]],outname=paste0("H",i,".kegg"),cutoff=0.05,heatmap = T)
}


avg_expr_B=apply(cpm(y),1,function(x) mean(x[which(group=="B")]))
avg_expr_H=list()
for (i in 1:10) {
  avg_expr_H[[i]]=apply(cpm(y),1,function(x) mean(x[which(group==paste0("H",i))]))
}

log_fc_H=matrix(NA,nrow = 10780,ncol = 10)

for (i in 1:10) {
  log_fc_H[,i]=log2(avg_expr_H[[i]]/avg_expr_B)
}

row.names(log_fc_H)=row.names(dat_B_H_filtered)

dat.d.m=c()
for ( i in 1:10){
  tmp.dat=log_fc_H[res.ID$flybase_gene_id,i]
  dat.d.m=cbind(dat.d.m,tmp.dat)
}

dat.d.m.pathview=dat.d.m
dat.d.m=as.data.frame(dat.d.m)
rownames(dat.d.m.pathview)=res.ID$entrezgene_id
dat.d.m.pathview=dat.d.m.pathview[!is.na(rownames(dat.d.m.pathview)),]
pathview(gene.data =dat.d.m.pathview, pathway.id = "dme00010", species = "dme",out.suffix="pathview")

kegg.gs.meta=kegg.gs[which(as.numeric(substr(names(kegg.gs),5,8))<1100)]

for (i in names(kegg.gs.meta)) {
  kegg.gs.meta[[i]]=kegg.gs.meta[[i]][kegg.gs.meta[[i]]%in%row.names(dat.d.m)]
}


kegg.gs.meta.heat.use=c()
for (i in names(kegg.gs.meta)) {
  temp=apply(dat.d.m[kegg.gs.meta[[i]],],2,mean)
  kegg.gs.meta.heat.use=rbind(kegg.gs.meta.heat.use,temp)
  }
row.names(kegg.gs.meta.heat.use)=names(kegg.gs.meta)
sort(apply(kegg.gs.meta.heat.use, 1, mean))

par(mfrow=c(3,3))
for (i in names(sort(apply(kegg.gs.meta.heat.use, 1, mean)))[1:5]) {
  temp=table(factor(table(unlist(lapply(sig_ID_dn,function(x) x[x%in%kegg.gs.meta[[i]]]))),levels=1:10))
  barplot(temp/sum(temp),ylim=c(0,1),main = i)
  text(9,0.8,length(kegg.gs.meta[[i]]))
}

colnames(kegg.gs.meta.heat.use)=paste0("rep",1:10)
pheatmap(kegg.gs.meta.heat.use,show_rownames = F,cluster_cols = F)

