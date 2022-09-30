# Rscript of correlation analysis and visualization
# author: Wei-Yun Lai
# date created: 2021.04.
# date last edited: 2022.08.19
# reviewed by Sheng-Kai Hsu 2022.08.21
# suggestion: need to rerun eveything as the binning method changes. Otherwise, we should change the causal analysis part.

rm(list = ls())
setwd('/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/project2/pleiotropy/')
####data input####
pleio_use_no01=read.table("./pleio_putatively_adaptive_gene.txt",sep = "\t",header = T,stringsAsFactors = F)
pleio_use_no01$f.value=-log(pleio_use_no01$f.value)


####replicate frequency spectrum (duplicated)####
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
sig_ID_all=list()

for (i in 1:10){
  my.contrasts_m=rep(0,length(unique(group_rep)))
  my.contrasts_m[c(1,i+1)]=c(-1,1)
  LRT_m=glmLRT(GLM_rep,contrast = my.contrasts_m)
  res_m=LRT_m$table
  row.names(res_m)=row.names(dat_B_H_filtered)
  res_m$padj=p.adjust(res_m$PValue,method = "BH")
  sig_ID_all[[i]]=row.names(res_m[res_m$padj<0.05,])
}


#png(filename = "/Users/weiyun/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/figure1b_v2.png")
barplot(table(table(unlist(sig_ID_all))),
        ylab = "Number of genes", xlab = "Number of populations",
        main = "Replicate frequency spectum")
#dev.off()

####parallelism####
#png("~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/figure1c.png",height = 10,width = 10,units = "cm",pointsize = 8,res = 600)
hist(pleio_use_no01$f.value,breaks = 20,col = "grey",
     main = "Distribution of evolutionary parallelism across genes",
     xlab = "evolutionary parallelism (log(1/F))")
#dev.off()

#png("~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/supplementaryfigure2.png",height = 10,width = 10,units = "cm",pointsize = 8,res = 600)
boxplot(pleio_use_no01$f.value~pleio_use_no01$rep_spe,col = "grey",ylim=c(-4.5,2),
     main = "",
     xlab = "number of populations",
     ylab = "evolutionary parallelism (log(1/F))")
#text(1:9,2,labels = c("c","ab","a","bc","c","d","d","de","e"))
#dev.off()

####index of pleiotropy/ancestral variance####

#grouping (tissue specificity)
ind=quantile(pleio_use_no01$ts_1.tau,seq(0,1,0.1))
group_tau=c()
for (i in 1:10) {
  group_tau[pleio_use_no01$ts_1.tau>=ind[i] & pleio_use_no01$ts_1.tau<ind[i+1]]=i
}

#grouping (connectivity)
ind=quantile(pleio_use_no01$connectivity,seq(0,1,0.1))
group_conn=c()
for (i in 1:10) {
  group_conn[pleio_use_no01$connectivity>=ind[i] & pleio_use_no01$connectivity<ind[i+1]]=i
}

#grouping (anc_v)
ind=quantile(pleio_use_no01$anc_v,seq(0,1,0.1))
group_anc=c()
for (i in 1:10) {
  group_anc[pleio_use_no01$anc_v>=ind[i] & pleio_use_no01$anc_v<ind[i+1]]=i
}


####correlation between connectivity and 1-tau####
cor(pleio_use_no01$ts_1.tau,pleio_use_no01$connectivity,method = "spearman")
#png(filename = "~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/sf1.png",height = 8.7,width = 8.7,res=600,pointsize = 6,units = "cm")
boxplot(pleio_use_no01$ts_1.tau~group_conn,
        xlab = "network connectivity",
        ylab = expression(1-tau))
text(2,0.8,labels = c("rho=0.52, p<0.001"),col="red")
#dev.off()

####correlation between anc_v and pleiotropy####
cor.test(log(pleio_use_no01$anc_v),pleio_use_no01$ts_1.tau,method = "spearman")
cor.test(log(pleio_use_no01$anc_v),pleio_use_no01$connectivity,method = "spearman")

#png(filename = "~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/figure4a.png",height = 13.05,width = 8.7,units = "cm",res = 600,pointsize = 8)
boxplot(log(as.numeric(paste(pleio_use_no01$anc_v)))~group_tau,xlab = "pleiotropy(1-tau)",ylab = "expression variation in the ancestral population",
        names=c("","20%","","40%","","60%","","80%","","100%"),col="grey")
text(8.5,0.25,labels = c("rho=-0.34, p<0.001"),col="red")
#dev.off()

#png(filename = "~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/sf3.png",height = 13.05,width = 8.7,units = "cm",res = 600,pointsize = 8)
boxplot(as.numeric(paste(pleio_use_no01$anc_v))~group_conn,
        names=c("20%","","40%","","60%","","80%","","100%"),xlab = "pleiotropy(connectivity)",log="y"
        ,ylab = "expression variance in ancestral population")
text(4,1.5,labels = c("rho=-0.37, p<0.001"),col="red")
#dev.off()

####correlation between anc_v and hetergeneity####
cor.test(log(pleio_use_no01$anc_v),pleio_use_no01$f.value,method = "spearman")

#png(filename = "~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/figure4b.png",height = 13.05,width = 8.7,units = "cm",res = 600,pointsize = 8)
boxplot(pleio_use_no01$f.value~group_anc,ylab = "evolutionary parallelism",xlab = "expression variation \nin the ancestral population",
        names=c("","20%","","40%","","60%","","80%","","100%"),col="grey",ylim=c(-4,2))
text(8.5,2,labels = c("rho=-0.24, p<0.001"),col="red")
#dev.off()


####correlation between heterogeneity and pleiotropy####
cor.test(pleio_use_no01$f.value,pleio_use_no01$ts_1.tau,method = "spearman")
cor.test(pleio_use_no01$f.value,pleio_use_no01$connectivity,method = "spearman")

#png(filename = "~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/figure2.png",height = 13.05,width = 8.7,units = "cm",res = 600,pointsize = 8)
boxplot(pleio_use_no01$f.value~group_tau,xlab = "pleiotropy(1-tau)",ylab = "evolutionary parallelism (log(1/F))",ylim=c(-4,2),
        names=c("","20%","","40%","","60%","","80%","","100%"),col="grey")
text(8.5,2,labels = "rho=0.21,p<0.001",col = "red")
#dev.off()

#boxplot(pleio_use_no01$ts_1.tau~pleio_use_no01$rep_spe,xlab = "number of replicates",ylab = "pleiotropy(1-tau)",
#        col="grey")

#png(filename = "~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_pleiotropy/sf2.png",height = 13.05,width = 8.7,units = "cm",res = 600,pointsize = 8)
boxplot(pleio_use_no01$f.value~group_conn,ylim=c(-4,2),
        names=c("20%","","40%","","60%","","80%","","100%"),xlab = "pleiotropy(connectivity)",
        ylab = "evolutionary parallelism",col="grey")
text(4,2,labels = "rho=0.10,p<0.001",col = "red")
#dev.off()


