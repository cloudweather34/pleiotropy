setwd("/Volumes/cluster/Wei-Yun/project2/GOdatabase/")

GO_res_table1=list()
for (i in 1:10){
  tmp=factor(as.integer(rownames(dat_B_H_filtered)%in%unlist(sig_ID_dn[[i]])))
  names(tmp)=rownames(dat_B_H_filtered)
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.classic",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  #  tmp_res$Fisher.classic=resTopGO.classic@score
  GO_res_table1[[i]]=tmp_res
}

GO_res_table1=lapply(GO_res_table1,function(x) {
  x$Fisher.classic[x$Fisher.classic=="< 1e-30"]=1e-30
  x$Fisher.classic.padj=p.adjust(x$Fisher.classic,method = "BH")
  return(x)})

GO_res_table2=list()
for (i in 1:10){
  tmp=factor(as.integer(rownames(dat_B_H_filtered)%in%unlist(sig_ID_up[[i]])))
  names(tmp)=rownames(dat_B_H_filtered)
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.classic",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  #  tmp_res$Fisher.classic=resTopGO.classic@score
  GO_res_table2[[i]]=tmp_res
}


GO_res_table2=lapply(GO_res_table2,function(x) {
  x$Fisher.classic[x$Fisher.classic=="< 1e-30"]=1e-30
  x$Fisher.classic.padj=p.adjust(x$Fisher.classic,method = "BH")
  return(x)})

ind=combinations(10,2,set=TRUE, repeats.allowed=FALSE)

GO.use.down=lapply(GO_res_table1,function(x) x$GO.ID[which(x$Fisher.classic.padj<0.05)])
GO.use.up=lapply(GO_res_table2,function(x) x$GO.ID[which(x$Fisher.classic.padj<0.05)])

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


sig_ID_dn_rs=list()
for (i in 1:10){
  sig_ID_dn_rs[[i]]=sample(rownames(res_m),sapply(sig_ID_dn,length)[i],replace = F)
}

GO_res_table3=list()
for (i in 1:10){
  tmp=factor(as.integer(rownames(dat_B_H_filtered)%in%unlist(sig_ID_dn_rs[[i]])))
  names(tmp)=rownames(dat_B_H_filtered)
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.classic",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  #tmp_res$Fisher.classic=resTopGO.classic@score
  GO_res_table3[[i]]=tmp_res
}

GO_res_table3=lapply(GO_res_table3,function(x) {
  x$Fisher.classic[x$Fisher.classic=="< 1e-30"]=1e-30
  x$Fisher.classic.padj=p.adjust(x$Fisher.classic,method = "BH")
  return(x)})

GO.use.down.rs=lapply(GO_res_table3,function(x) x$GO.ID[which(x$Fisher.classic.padj<0.05)])

ja_go_dn_rs=c()
for (i in 1:45) {
  inte=length(intersect(unlist(GO.use.down.rs[ind[i,1]]),unlist(GO.use.down.rs[ind[i,2]])))
  uni=length(union(unlist(GO.use.down.rs[ind[i,1]]),unlist(GO.use.down.rs[ind[i,2]])))
  ja_go_dn_rs[i]=inte/uni
}

ja_dn_rs=c()
for (i in 1:45) {
  inte=length(intersect(unlist(sig_ID_dn_rs[ind[i,1]]),unlist(sig_ID_dn_rs[ind[i,2]])))
  uni=length(union(unlist(sig_ID_dn_rs[ind[i,1]]),unlist(sig_ID_dn_rs[ind[i,2]])))
  ja_dn_rs[i]=inte/uni
}

boxplot(ja_dn_rs,ja_go_dn_rs)
####TopGO####
temp=as.factor(dat_B_H$X%in%dat_B_H$X[1000:1100])
names(temp)=dat_B_H$X
tgd1=new( "topGOdata", ontology="BP", allGenes = temp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
graph(tgd1)
buildLevels(graph(tgd1))
GO_list=as.list(buildLevels(graph(tgd1))$level2nodes)

####GOsim####
GO_all_gene_ID_length=sapply(genesInTerm(object = tgd1),length)


consi=strsplit2("GO:0006091,GO:0006163,GO:0006164,GO:0006754,GO:0009144,GO:0009145,GO:0009150,GO:0009152,GO:0009199,GO:0009201,GO:0009205,GO:0009206,GO:0009259,GO:0009260,GO:0010927,GO:0015672,GO:0015985,GO:0015986,GO:0019693,GO:0030239,GO:0031032,GO:0034220,GO:0040003,GO:0042335,GO:0042692,GO:0045214,GO:0046034,GO:0046390,GO:0051146,GO:0055001,GO:0055002,GO:0072521,GO:0098655,GO:0098660,GO:0098662",",")[1,]


lev=19-sapply(consi,function(x) which(sapply(GO_list,function(y) any(grepl(x,y)))))


jac=as.numeric(strsplit2("0.62274606972906,0.621044848154709,0.719426899484108,0.816084656084656,0.816084656084656,0.816084656084656,0.621044848154709,0.719426899484108,0.816084656084656,0.816084656084656,0.816084656084656,0.816084656084656,0.618045537330518,0.712194749059738,0.570726127806413,0.538792584945673,0.802543752543753,0.802543752543753,0.618045537330518,0.727325534942045,0.707463559292998,0.477521330761709,0.522974854473171,0.516691584092846,0.634360647242186,0.678156070509012,0.618985224898413,0.712194749059738,0.650234849560173,0.673262750694481,0.673262750694481,0.623726306965525,0.502700864858975,0.525261708701279,0.512994645113815"
,",")[1,])
