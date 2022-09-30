rm(list=ls())
library(MetaboAnalystR)
library(globaltest)
library(edgeR)
library(metaseqR)
setwd("/Volumes/Temp2/shengkai/WY/")

####transcriptome: globaltest####
#BiocManager::install("globaltest")
#library(golubEsets)
#data(Golub_Train)
#library(vsn)
#exprs(Golub_Train) <- exprs(vsn2(Golub_Train))

dat_B_H=read.csv(file = "/Volumes/Temp2/shengkai/WY/monster_fullset_readcounts.csv",sep = ";")
colnames(dat_B_H)[1]="name"
dat_B_H_filtered=dat_B_H[apply(cpm(dat_B_H[,-1]),1,function(x){!any(x<=1)}),]

row.names(dat_B_H_filtered)=dat_B_H_filtered[,1]
dat_B_H_filtered=dat_B_H_filtered[,-1]
group_rep=paste(substr(colnames(dat_B_H_filtered),1,1),substr(colnames(dat_B_H_filtered),8,8),sep = "_")
group_rep[c(11,12,23,33,34)]="B"

ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
conv_ID=ID_converter(ID = rownames(dat_B_H_filtered),db = ensembl,attributes = c("flybase_gene_id","entrezgene_id"),filters = "flybase_gene_id")

dat_use=dat_B_H_filtered[na.omit(conv_ID)$flybase_gene_id,]
rownames(dat_use)=na.omit(conv_ID)$entrezgene_id
y=DGEList(counts=dat_use,group = substr(colnames(dat_B_H_filtered),1,1))
y=calcNormFactors(y)

glob.test.res=list()
for(i in 1:10){
  idx=group_rep%in%c("B",paste0("H_",i-1))
  meta.table=data.frame(group = as.factor(group_rep[idx]),
                        row.names = rownames(y$samples)[idx])
  annot=data.frame(labelDescription=c("Factor levels"))
  annot_factors=AnnotatedDataFrame(meta.table,annot)
  expr.set=ExpressionSet(log(cpm(y)[,idx]),annot_factors)
  gt.options(transpose=TRUE)
  res=gt(ALL.AML, Golub_Train)
  gt(group, expr.set)
  glob.test.res[[i]]=gtKEGG(group, expr.set,annotation = "org.Dm.eg.db",multtest = "BH")
}
lapply(glob.test.res,function(x) dim(x@result))

sig_pathway_ge=lapply(glob.test.res,function(x) x@extra$alias[x@extra$'BH'<0.05])

ind=combinations(10,2,set=TRUE, repeats.allowed=FALSE)

ja_ge_path=c()
for (j in 1:45) {
  inte=length(intersect(unlist(sig_pathway_ge[ind[j,1]]),unlist(sig_pathway_ge[ind[j,2]])))
  uni=length(union(unlist(sig_pathway_ge[ind[j,1]]),unlist(sig_pathway_ge[ind[j,2]])))
  ja_ge_path[j]=inte/uni
}


####metabolome: metaboanalystR: QSEA####
dat=read.delim("./met_raw_sk_modified.txt",sep = "\t",stringsAsFactors = F,header = F)
dat_rmdupl=dat[!duplicated(dat$V1),]
dat_rmdupl_out=t(dat_rmdupl)
#dat_rmdupl_log[-1,-1:-2]=log(as.numeric(dat_rmdupl_log[-1,-1:-2]))
dat_rmdupl_out[-1,-1:-2]=as.numeric(dat_rmdupl_out[-1,-1:-2])
dat_rmdupl_out[1,][51]="C14530"
for (i in 1:6){
  print(dat_rmdupl_out[c(1:6,6+(i*3-2):(i*3)),1])
  tmp_dat=dat_rmdupl_out[c(1:6,6+(i*3-2):(i*3)),]
  write.table(tmp_dat,paste0("./met_data_metaboanalystR",i+3,".txt"),quote = F,col.names = F,row.names = F,sep = "\t")
}

write.table(dat_rmdupl_out,"./met_data_metaboanalystR.txt",quote = F,col.names = F,row.names = F,sep = "\t")

res_tab_qsea=list()
for (i in 1:6){
  rm(mSet)
  mSet=InitDataObjects("conc", "msetqea", FALSE)
  mSet=Read.TextData(mSet, paste0("./met_data_metaboanalystR",i+3,".txt"), "rowu", "disc");
  mSet=CrossReferencing(mSet, "kegg")
  mSet=CreateMappingResultTable(mSet)
  mSet=SanityCheckData(mSet)
  mSet=ReplaceMin(mSet)
  mSet=PreparePrenormData(mSet)
  mSet=Normalization(mSet, "NULL", "LogNorm", "NULL", "B1", ratio=FALSE, ratioNum=20)
  mSet=PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet=PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  mSet=SetMetabolomeFilter(mSet, F)
  mSet=SetCurrentMsetLib(mSet, "kegg_pathway", 1)
  mSet=CalculateGlobalTestScore(mSet)
  mSet=PlotQEA.Overview(mSet, "qea_0_", "bar", "png", 72, width=NA)
  res_tab_qsea[[i]]=mSet$analSet$qea.mat
}

sig_pathway=lapply(res_tab_qsea,function(x) rownames(x)[x$'Raw p'<0.05])

library(gtools)
ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE)

ja_met_path=c()
for (j in 1:15) {
  inte=length(intersect(unlist(sig_pathway[ind[j,1]]),unlist(sig_pathway[ind[j,2]])))
  uni=length(union(unlist(sig_pathway[ind[j,1]]),unlist(sig_pathway[ind[j,2]])))
  ja_met_path[j]=inte/uni
}
####intergrated analysis####
res_tab_comb=list()
for( i in 1:6){
  tmp=data.frame(pathway=glob.test.res[[i+4]]@extra$alias,p_gene=glob.test.res[[i+4]]@result[,1])
  for (j in 1:dim(tmp)[1]){
    if(tmp$pathway[j]%in%rownames(res_tab_qsea[[i]])) {
      tmp$p_meta[j]=res_tab_qsea[[i]][tmp$pathway[j],]$`Raw p`
      }
    else tmp$p_meta[j]=1
  }
  tmp$p_comb=fisher.method(pvals = cbind(tmp$p_gene,tmp$p_meta),p.corr = "BH")[,4]
  tmp$FDR_comb=p.adjust(tmp$p_comb,method = "BH")
  res_tab_comb[[i]]=tmp[order(tmp$FDR_comb),]                 
}

res_tab_comb=lapply(res_tab_comb,function(x) {
  x$kegg_map=rownames(x)
  return(x)})

for (i in 1:6){
  write.table(res_tab_comb[[i]],paste0("./joint_pathway_analysis_R",i+3,".txt"),quote = F,
              row.names = F,col.names = T,sep="\t")
}

sig_path_comb=lapply(res_tab_comb,function(x) x$pathway[x$FDR_comb<0.05])

ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE)

ja_comb_path=c()
for (j in 1:15) {
  inte=length(intersect(unlist(sig_path_comb[ind[j,1]]),unlist(sig_path_comb[ind[j,2]])))
  uni=length(union(unlist(sig_path_comb[ind[j,1]]),unlist(sig_path_comb[ind[j,2]])))
  ja_comb_path[j]=inte/uni
}



