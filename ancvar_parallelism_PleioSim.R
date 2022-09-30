# Rscript to process the simulated data for the relationship between ancestral variance and evolutionary parallelism
# date created: 2022.04.10
# date last edited: 2022.08.21

rm(list = ls())
library(poolSeq)
setwd("/Volumes/Temp3/Wei-Yun/Wei-Yun/simulation/")

####read sync and effect sizes####
filename.100=list.files("./pleiotropyS2_100/",pattern = ".gz")[1:99] #toy dataset
es.name.100=list.files("./pleiotropyS2_100/",pattern = "effectsize")[1:99] #toy dataset

sync.list.100=lapply(filename.100,function(x) {
  print(x)
  read.sync(paste0("./pleiotropyS2_100/",x),
            gen=rep(c(0,1,5,10,15,20,25,35,50,75,100,200),10),
            repl=rep(1:10,each=12),polarization = "minor")})

es.list.100=lapply(es.name.100,function(x)
  read.table(paste0("./pleiotropyS2_100/",x),
             header = F,sep="\t",stringsAsFactors = F))

####calculate allele frequency and changes####
af.list.100=lapply(sync.list.100, function(x){
  af(x,gen=rep(c(0,1,5,10,15,20,25,35,50,75,100,200),10),repl=rep(1:10,each=12))
})

#creating 3 sub-replicates by random sampling 50 individuals
af.list.100.s1=lapply(af.list.100,function(x) apply(x,c(1,2),function(a) sum(sample(c(0,1),50,prob = c(1-a,a),replace = T))/50))
af.list.100.s2=lapply(af.list.100,function(x) apply(x,c(1,2),function(a) sum(sample(c(0,1),50,prob = c(1-a,a),replace = T))/50))
af.list.100.s3=lapply(af.list.100,function(x) apply(x,c(1,2),function(a) sum(sample(c(0,1),50,prob = c(1-a,a),replace = T))/50))

sf.list.100=lapply(af.list.100, function(x){
  x[,1]
}) # starting frequency
afc.list.100=lapply(af.list.100, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
}) # AFC by generations

afc.list.100.s1=lapply(af.list.100.s1, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
}) #AFC in repl.1
afc.list.100.s2=lapply(af.list.100.s2, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
}) #AFC in repl.2
afc.list.100.s3=lapply(af.list.100.s3, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
}) #AFC in repl.3
afc.list.100.sample=list(afc.list.100.s1,afc.list.100.s2,afc.list.100.s3) #combine into list

####process effect size file and rescale to expression####
es.processed.list.100=lapply(es.list.100,function(x){
  tmp=x[,4]
  names(tmp)=paste(x[,1],x[,2],sep=".")
  tmp[1:5]=tmp[1:5]/20
  tmp[6:20]=tmp[6:20]*15/100
  tmp[21:50]=tmp[21:50]*30/100
  tmp[51:100]=tmp[51:100]*50/100
  return(tmp)
})

#### calculate ancestral genetic variance and difference in genetic value (expression)####
#four genes were calculated indepently
anc.v.1=c()
anc.v.2=c()
anc.v.3=c()
anc.v.4=c()
d.gv.1=c()
d.gv.2=c()
d.gv.3=c()
d.gv.4=c()
for(i in 1:length(afc.list.100)){
  # assign the contributing loci to each gene
  idx=rep(0,100)
  idx[names(sf.list.100[[i]])%in%names(es.processed.list.100[[i]])[1:5]]=1
  idx[names(sf.list.100[[i]])%in%names(es.processed.list.100[[i]])[6:20]]=2
  idx[names(sf.list.100[[i]])%in%names(es.processed.list.100[[i]])[21:50]]=3
  idx[names(sf.list.100[[i]])%in%names(es.processed.list.100[[i]])[51:100]]=4
  
  #anc. variance = 2pqa^2
  anc.v.1=c(anc.v.1,sum(2*1e4*sf.list.100[[i]][idx==1]*(1-sf.list.100[[i]][idx==1])*es.processed.list.100[[i]][names(sf.list.100[[i]])][idx==1]^2))
  anc.v.2=c(anc.v.2,sum(2*1e4*sf.list.100[[i]][idx==2]*(1-sf.list.100[[i]][idx==2])*es.processed.list.100[[i]][names(sf.list.100[[i]])][idx==2]^2))
  anc.v.3=c(anc.v.3,sum(2*1e4*sf.list.100[[i]][idx==3]*(1-sf.list.100[[i]][idx==3])*es.processed.list.100[[i]][names(sf.list.100[[i]])][idx==3]^2))
  anc.v.4=c(anc.v.4,sum(2*1e4*sf.list.100[[i]][idx==4]*(1-sf.list.100[[i]][idx==4])*es.processed.list.100[[i]][names(sf.list.100[[i]])][idx==4]^2))

  # difference in genetic value = delta(p) * a
  d.tmp1=c()
  d.tmp2=c()
  d.tmp3=c()
  d.tmp4=c()
  for (j in 1:3){
    tmp1=apply(afc.list.100.sample[[j]][[i]],2,function(x) {
      sum(x[idx==1]*100*es.processed.list.100[[i]][rownames(afc.list.100[[i]])][idx==1])
    })
    d.tmp1=c(d.tmp1,tmp1)
    tmp2=apply(afc.list.100.sample[[j]][[i]],2,function(x) {
      sum(x[idx==2]*100*es.processed.list.100[[i]][rownames(afc.list.100[[i]])][idx==2])
    })
    d.tmp2=c(d.tmp2,tmp2)
    tmp3=apply(afc.list.100.sample[[j]][[i]],2,function(x) {
      sum(x[idx==3]*100*es.processed.list.100[[i]][rownames(afc.list.100[[i]])][idx==3])
    })
    d.tmp3=c(d.tmp3,tmp3)
    tmp4=apply(afc.list.100.sample[[j]][[i]],2,function(x) {
      sum(x[idx==4]*100*es.processed.list.100[[i]][rownames(afc.list.100[[i]])][idx==4])
    })
    d.tmp4=c(d.tmp4,tmp4)
  }
  d.gv.1=rbind(d.gv.1,d.tmp1)
  d.gv.2=rbind(d.gv.2,d.tmp2)
  d.gv.3=rbind(d.gv.3,d.tmp3)
  d.gv.4=rbind(d.gv.4,d.tmp4)  
}

f.dgv.1=apply(d.gv.1,1,fstat)
f.dgv.2=apply(d.gv.2,1,fstat)
f.dgv.3=apply(d.gv.3,1,fstat)
f.dgv.4=apply(d.gv.4,1,fstat)

par(mfrow=c(1,2))
boxplot(anc.v.1,anc.v.2,anc.v.3,anc.v.4,xlab="# of loci",ylab="Variance of ancestral expression",
        names=c(5,15,30,50))
boxplot(log(1/f.dgv.1),log(1/f.dgv.2),log(1/f.dgv.3),log(1/f.dgv.4),xlab="# of loci",ylab="Variance of expression changes",
        names=c(5,15,30,50))

cor.test(c(anc.v.1,anc.v.2,anc.v.3,anc.v.4),log(1/c(f.dgv.1,f.dgv.2,f.dgv.3,f.dgv.4)),method = "spearman")
fm=lm(log(1/c(f.dgv.1,f.dgv.2,f.dgv.3,f.dgv.4))~c(anc.v.1,anc.v.2,anc.v.3,anc.v.4))
png("~/Dropbox (PopGen)/backup/Wei-Yun/project2/pleiotropy/ancVar2Parallelism_simulation.png",units = "cm",
    width = 8.7,height = 8.7,res=600,pointsize = 8)
par(mar=c(5,5,2,2))
plot(c(anc.v.1,anc.v.2,anc.v.3,anc.v.4),log(1/c(f.dgv.1,f.dgv.2,f.dgv.3,f.dgv.4)),
     ylab="Gene expression evolution parallelism",xlab="Variance of ancestral expression")
abline(fm$coefficients,col="red")
text(12,4, labels=expression(rho==-0.26),col="red")
text(12,3.5,labels = expression(p==1.77e-7),col="red")
dev.off()

gpf.name=list.files("./opt1_0.5ef_10rep_100/",pattern = ".gpf")[1:100]
gpf.list=lapply(gpf.name,function(x) read.table(paste0("./opt1_0.5ef_10rep_100/",x),header = F))

avg_gv=sapply(gpf.list,function(x) tapply(x[,4],paste(x[,1],x[,2]),mean))


####from mimhap####
numericGeno=function(x){
  out=x[,-c(1:4)]
  allele=limma::strsplit2(x[,4],split ="/")
  out=apply(out,2,function(a) {
    res=rep(1,length(a))
    res[a==paste0(allele[,2],allele[,2])]=0
    res[a==paste0(allele[,1],allele[,1])]=2
    return(res)
  })
  colnames(out)=1:300
  rownames(out)=paste(x[,1],x[,2],sep = ".")
  return(out)
}
fstat=function(x){
  fm=lm(x~gl(10,k = 1,length = 30))
  stat=anova(fm)$'F value'[1]
  return(stat)
}
es.name.100=list.files("./pleiotropyS1_100/",pattern = "effectsize",full.names = T)[1:99] #toy dataset
es.list.100=lapply(es.name.100,function(x)
  read.table(x,
             header = F,sep="\t",stringsAsFactors = F))
es.processed.list.100=lapply(es.list.100,function(x){
  tmp=x[,4]
  names(tmp)=paste(x[,1],x[,2],sep=".")
  tmp[1:5]=tmp[1:5]/20
  tmp[6:20]=tmp[6:20]*15/100
  tmp[21:50]=tmp[21:50]*30/100
  tmp[51:100]=tmp[51:100]*50/100
  return(tmp)
})

hap.dir=list.dirs("./pleiotropyS1_100/")[2:100]
g1.hap.name=lapply(hap.dir,function(x) list.files(x,pattern = "g1.mimhap.gz",full.names = T))
g1.hap.list=lapply(g1.hap.name,function(x) lapply(x,function(a) read.table(a,header = F)))
g1.num.list=lapply(g1.hap.list,function(x) lapply(x,numericGeno))

g100.hap.name=lapply(hap.dir,function(x) list.files(x,pattern = "g100.mimhap.gz",full.names = T))
g100.hap.list=lapply(g100.hap.name,function(x) lapply(x,function(a) read.table(a,header = F)))
g100.num.list=lapply(g100.hap.list,function(x) lapply(x,numericGeno))

set.seed(100)
idxLoci=rep(1:4,c(5,15,30,50))

g1.gv.list=list()
g100.gv.list=list()
for(i in 1:99){
  idxLoci=rep(0,100)
  idxLoci[rownames(g1.num.list[[i]])%in%names(es.processed.list.100[[i]])[1:5]]=1
  idxLoci[rownames(g1.num.list[[i]])%in%names(es.processed.list.100[[i]])[6:20]]=2
  idxLoci[rownames(g1.num.list[[i]])%in%names(es.processed.list.100[[i]])[21:50]]=3
  idxLoci[rownames(g1.num.list[[i]])%in%names(es.processed.list.100[[i]])[51:100]]=4
  table(idxLoci)
  g1.gv.list[[i]]=lapply(g1.num.list[[i]],function(a) apply(a,2,function(b) b*es.processed.list.100[[i]][rownames(a)]))
  g100.gv.list[[i]]=lapply(g100.num.list[[i]],function(a) apply(a,2,function(b) b*es.processed.list.100[[i]][rownames(a)]))
}

g1.trgv.list=lapply(g1.gv.list,function(x) lapply(x,function(a) apply(a,2,function(b) tapply(b,idxLoci,sum))))
g100.trgv.list=lapply(g100.gv.list,function(x) lapply(x,function(a) apply(a,2,function(b) tapply(b,idxLoci,sum))))

g1.trgv.list.s1=lapply(g1.trgv.list,function(x) sapply(x,function(a) rowSums(a[,sample(1:300,50)])/50))
g1.trgv.list.s2=lapply(g1.trgv.list,function(x) sapply(x,function(a) rowSums(a[,sample(1:300,50)])/50))
g1.trgv.list.s3=lapply(g1.trgv.list,function(x) sapply(x,function(a) rowSums(a[,sample(1:300,50)])/50))
g100.trgv.list.s1=lapply(g100.trgv.list,function(x) sapply(x,function(a) rowSums(a[,sample(1:300,50)])/50))
g100.trgv.list.s2=lapply(g100.trgv.list,function(x) sapply(x,function(a) rowSums(a[,sample(1:300,50)])/50))
g100.trgv.list.s3=lapply(g100.trgv.list,function(x) sapply(x,function(a) rowSums(a[,sample(1:300,50)])/50))

g1.trgv.list.ss=list()
g100.trgv.list.ss=list()
for(i in 1:99){
  g1.trgv.list.ss[[i]]=cbind(g1.trgv.list.s1[[i]],g1.trgv.list.s2[[i]],g1.trgv.list.s3[[i]])
  g100.trgv.list.ss[[i]]=cbind(g100.trgv.list.s1[[i]],g100.trgv.list.s2[[i]],g100.trgv.list.s3[[i]])
}
g1.trgv.array.ss=simplify2array(g1.trgv.list.ss)
g100.trgv.array.ss=simplify2array(g100.trgv.list.ss)

p.val=array(dim = c(4,10,99))
for (i in 1:99){
  for (j in 1:10){
    for (k in 1:4){
      p.val[k,j,i]=t.test(g1.trgv.array.ss[k,c(j,j+10,j+20),i],g100.trgv.array.ss[k,c(j,j+10,j+20),i])$p.value
    }
  }
}
idx=apply(p.val<0.05/99,3,function(x) rowSums(x)>1)

d.trgv=g100.trgv.array.ss - g1.trgv.array.ss[,rep(1:10,3),]
f.d.trgv=apply(d.trgv,c(1,3),fstat)
boxplot(t(f.d.trgv))

####alternative analysis
filename.100=list.files("./opt1_2.5ef_10rep_100/",pattern = ".gz")[1:100] #toy dataset
filename.30=list.files("./opt1_2.5ef_10rep_10/",pattern = ".gz")[1:100] #toy dataset

es.name.100=list.files("./opt1_2.5ef_10rep_100/",pattern = "effectsize")[1:100] #toy dataset
es.name.30=list.files("./opt1_2.5ef_10rep_10/",pattern = "effectsize")[1:100] #toy dataset

sync.list.100=lapply(filename.100,function(x) 
  read.sync(paste0("./opt1_2.5ef_10rep_100/",x),
            gen=rep(c(0,1,5,10,15,20,25,35,50,75,100,200),10),
            repl=rep(1:10,each=12),polarization = "minor"))
sync.list.30=lapply(filename.30,function(x) 
  read.sync(paste0("./opt1_2.5ef_10rep_10/",x),
            gen=rep(c(0,1,5,10,15,20,25,35,50,75,100,200),10),
            repl=rep(1:10,each=12),polarization = "minor"))

es.list.100=lapply(es.name.100,function(x)
  read.table(paste0("./opt1_2.5ef_10rep_100/",x),
             header = F,sep="\t",stringsAsFactors = F))
es.list.30=lapply(es.name.30,function(x)
  read.table(paste0("./opt1_2.5ef_10rep_10/",x),
             header = F,sep="\t",stringsAsFactors = F))

af.list.100=lapply(sync.list.100, function(x){
  af(x,gen=rep(c(0,1,5,10,15,20,25,35,50,75,100,200),10),repl=rep(1:10,each=12))
})
af.list.30=lapply(sync.list.30, function(x){
  af(x,gen=rep(c(0,1,5,10,15,20,25,35,50,75,100,200),10),repl=rep(1:10,each=12))
})

af.list.100.sample=list()
i=1
while (i<=3){
  af.list.100.sample[[i]]=lapply(af.list.100,function(x) apply(x,c(1,2),function(a) sum(sample(c(0,1),100,prob = c(1-a,a),replace = T))/100))
  i=i+1
}
af.list.30.sample=list()
i=1
while (i<=3){
  af.list.30.sample[[i]]=lapply(af.list.30,function(x) apply(x,c(1,2),function(a) sum(sample(c(0,1),100,prob = c(1-a,a),replace = T))/100))
  i=i+1
}

sf.list.100=lapply(af.list.100, function(x){
  x[,1]
})
sf.list.30=lapply(af.list.30, function(x){
  x[,1]
})
afc.list.100=lapply(af.list.100, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
})
afc.list.30=lapply(af.list.30, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
})

afc.list.100.sample=lapply(af.list.100.sample,function(a) lapply(a, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
}))

afc.list.30.sample=lapply(af.list.30.sample,function(a) lapply(a, function(x){
  tmp=x[,c(1,seq(11,120,12))]
  afc=tmp-tmp[,1]
  return(afc[,-1])
}))

es.processed.list.100=lapply(es.list.100,function(x){
  tmp=x[,4]
  names(tmp)=paste(x[,1],x[,2],sep=".")
  return(tmp)
})

es.processed.list.30=lapply(es.list.30,function(x){
  tmp=x[,4]
  names(tmp)=paste(x[,1],x[,2],sep=".")
  return(tmp)
})

anc.v.100=c()
anc.v.30=c()
d.gv.all.100=c()
d.gv.all.30=c()
d.gv.100=c()
d.gv.30=c()
for(i in 1:length(afc.list.100)){
  anc.v.100=c(anc.v.100,sum(2*sf.list.100[[i]]*(1-sf.list.100[[i]])*es.processed.list.100[[i]]^2))
  anc.v.30=c(anc.v.30,sum(2*sf.list.30[[i]]*(1-sf.list.30[[i]])*es.processed.list.30[[i]]^2))
  tmp.100=apply(afc.list.100[[i]],2,function(x) {
    sum(x*es.processed.list.100[[i]][rownames(afc.list.100[[i]])])
  })
  d.gv.all.100=rbind(d.gv.all.100,tmp.100)
  tmp.30=apply(afc.list.30[[i]],2,function(x) {
    sum(x*es.processed.list.30[[i]][rownames(afc.list.30[[i]])])
  })
  d.gv.all.30=rbind(d.gv.all.30,tmp.30)
  
  d.tmp1=c()
  d.tmp2=c()
  for (j in 1:3){
    tmp1=apply(afc.list.100.sample[[j]][[i]],2,function(x) {
      sum(x*es.processed.list.100[[i]][rownames(afc.list.100[[i]])])
    })
    d.tmp1=c(d.tmp1,tmp1)
    tmp2=apply(afc.list.30.sample[[j]][[i]],2,function(x) {
      sum(x*es.processed.list.30[[i]][rownames(afc.list.30[[i]])])
    })
    d.tmp2=c(d.tmp2,tmp2)
      }
  d.gv.100=rbind(d.gv.100,d.tmp1)
  d.gv.30=rbind(d.gv.30,d.tmp2)
}

v.dgv.100=apply(d.gv.100,1,fstat)
v.dgv.30=apply(d.gv.30,1,fstat)

par(mfrow=c(1,2))
boxplot(anc.v.30,anc.v.100,xlab="# of loci",ylab="Variance of ancestral expression",
        names=c(30,100))
boxplot(v.dgv.30,v.dgv.100,xlab="# of loci",ylab="Variance of expression changes",
        names=c(30,100))

plot(c(anc.v.30,anc.v.100),c(v.dgv.30,v.dgv.100),
     ylab="Variance of expression changes",xlab="Variance of ancestral expression")


adapt_cv_0.5ef_5=c()
adapt_sd_0.5ef_5=c()
adapt_r_0.5ef_5=c()
for (i in 1:length(filename.0.5ef)){
  snp.idx=paste(es.list.0.5ef[[i]][,1],es.list.0.5ef[[i]][,2],sep = ".")[abs(es.list.0.5ef[[i]]$V4)>1/5/2]
  af_all=af(sync.list.0.5ef[[i]],repl = 1:10,gen = c(1,5,10,15,20,25,35,50,75,100,200))
  gen.idx=c(grep("F1.R",colnames(af_all)),grep("F200",colnames(af_all)))
  AFC_matrix=t(apply(af_all[,gen.idx],1,function(x) tapply(x,rep(1:10,2),function(a) a[2]-a[1])))
  snp.idx2=names(which(apply(AFC_matrix,1,function(x) any(x>0.1))))
  snp.idx3=intersect(snp.idx,snp.idx2)
  print(length(snp.idx3))
  AFC_matrix=AFC_matrix[snp.idx3,]
  if (length(snp.idx3)>1){
    cv_AFC=apply(AFC_matrix,1,function(a) sd(a)/abs(mean(a)))
    sd_AFC=apply(AFC_matrix,1,sd)
    r_AFC=apply(af_all[snp.idx3,],1,function(x) mean(cor(matrix(x,11,10,byrow = F)),na.rm = T))
  }
  else{
    cv_AFC=sd(AFC_matrix)/abs(mean(AFC_matrix))
    sd_AFC=sd(AFC_matrix)
    r_AFC=mean(cor(matrix(af_all[snp.idx3,],11,10,byrow = F)),na.rm = T)
  }
  adapt_cv_0.5ef_5=c(adapt_cv_0.5ef_5,mean(cv_AFC))
  adapt_sd_0.5ef_5=c(adapt_sd_0.5ef_5,mean(sd_AFC))
  adapt_r_0.5ef_5=c(adapt_r_0.5ef_5,mean(r_AFC))
}

adapt_cv_100ef_5=c()
adapt_sd_100ef_5=c()
adapt_r_100ef_5=c()
for (i in 1:length(filename.100ef)){
  snp.idx=paste(es.list.100ef[[i]][,1],es.list.100ef[[i]][,2],sep = ".")[abs(es.list.100ef[[i]]$V4)>1/5/2]
  af_all=af(sync.list.100ef[[i]],repl = 1:10,gen = c(1,5,10,15,20,25,35,50,75,100,200))
  gen.idx=c(grep("F1.R",colnames(af_all)),grep("F200",colnames(af_all)))
  AFC_matrix=t(apply(af_all[,gen.idx],1,function(x) tapply(x,rep(1:10,2),function(a) a[2]-a[1])))
  snp.idx2=names(which(apply(AFC_matrix,1,function(x) any(x>0.1))))
  snp.idx3=intersect(snp.idx,snp.idx2)
  print(length(snp.idx3))
  AFC_matrix=AFC_matrix[snp.idx3,]
  if (length(snp.idx3)>1){
    cv_AFC=apply(AFC_matrix,1,function(a) sd(a)/abs(mean(a)))
    sd_AFC=apply(AFC_matrix,1,sd)
    r_AFC=apply(af_all[snp.idx3,],1,function(x) mean(cor(matrix(x,11,10,byrow = F)),na.rm = T))
  }
  else{
    cv_AFC=sd(AFC_matrix)/abs(mean(AFC_matrix),na.rm = T)
    sd_AFC=sd(AFC_matrix)
    r_AFC=mean(cor(matrix(af_all[snp.idx3,],11,10,byrow = F)))
  }
  adapt_cv_100ef_5=c(adapt_cv_100ef_5,mean(cv_AFC))
  adapt_sd_100ef_5=c(adapt_sd_100ef_5,mean(sd_AFC))
  adapt_r_100ef_5=c(adapt_r_100ef_5,mean(r_AFC))
}


par(mfrow=c(1,1))
#boxplot(adapt_sd_0.5ef_5,adapt_sd_100ef_5)
#boxplot(adapt_cv_0.5ef_5,adapt_cv_100ef_5)
#boxplot(adapt_r_0.5ef_5,adapt_r_100ef_5,ylab="Pearson's r of AFC",names())
#wilcox.test(adapt_sd_0.5ef_5,adapt_sd_100ef_5)
#wilcox.test(adapt_cv_0.5ef_5,adapt_cv_100ef_5)
#wilcox.test(adapt_r_0.5ef_5,adapt_r_100ef_5)

png(filename = "~/Dropbox (PopGen)/application_for_popgen_vienna/manuscript/figure2_S7.png",width = 8.7,height = 8.7,pointsize = 10,units = "cm",res = 600)
boxplot(adapt_r_0.5ef_5,adapt_r_100ef_5,
        xaxt="n",ylab="Pearson's r of AFC")
axis(1,at=1:2,padj=0.5,
     labels = paste(c("alpha=0.5","alpha=100"),
                    "5 loci",sep = "\n"))
dev.off()

