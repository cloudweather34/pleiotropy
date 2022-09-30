setwd("/Volumes/cluster/Wei-Yun/project2/")


net.con=read.table("flynet_supervised_0.6.txt",header=F,stringsAsFactors = F)


ind=unique(net.con$V1)
out_degree=as.data.frame(matrix(NA,617,2))
out_degree[,1]=as.character(unique(net.con$V1))
for (i in 1:617) {
  out_degree[i,2]=sum(net.con$V1==ind[i])
}
colnames(out_degree)=c("FBgn","out_de")


hist(as.numeric(paste(out_degree[,2])))

ind=unique(net.con$V2)
in_degree=as.data.frame(matrix(NA,12286,2))
in_degree[,1]=as.character(unique(net.con$V2))
for (i in 1:12286) {
  in_degree[i,2]=sum(net.con$V2==ind[i])
}
colnames(in_degree)=c("FBgn","in_de")

hist(as.numeric(paste(in_degree[,2])))

pleio_bg_in_tau=merge(aov_tau,in_degree,by = "FBgn",all = T)
pleio_evo=merge(pleio_bg_in_tau,out_degree,by="FBgn",all = T)

class(tau_bg$FBgn)
class(in_degree$FBgn)
class(out_degree$FBgn)

pleio_evo$conn=apply(pleio_evo[,c(7,8)],1,function(x) sum(x, na.rm = T))
plot(as.numeric(paste(pleio_bg$F_value)),pleio_evo$conn,log="y",xlim=c(0,50))
plot(pleio_evo$tau_male_4,pleio_evo$conn,log="y")
cor(pleio_evo$tau_male_4,pleio_evo$conn,use = "complete.obs",method = "spearman")

hist(log(pleio_evo$conn),xlab = "ln(connectivity)",main="distribtion of connectivity")

ind=seq(0,1,0.1)
group_tau=c()
for (i in 1:10) {
  group_tau[pleio_evo$tau_male_4>ind[i] & pleio_evo$tau_male_4<ind[i+1]]=i
}

boxplot(log(pleio_evo$conn)~group_tau,names=seq(0.15,0.95,0.1)
        ,xlab = "tissue specificity (group median)",ylab = "ln(connectivity)")


ind_conn=seq(0,10,2)
group_conn=c()
for (i in 1:5) {
  group_conn[log(pleio_evo$conn)>ind_conn[i] & log(pleio_evo$conn)<ind_conn[i+1]]=i
}


boxplot(log(as.numeric(paste(pleio_evo$F_value)))~group_conn,ylim=c(0,5),
        names=seq(1,9,2),xlab = "ln(connectivity) (group median)",
        ylab = "ln(F-value)")



ind=seq(0,1,0.2)
group_tau=c()
for (i in 1:5) {
  group_tau[pleio_evo$tau_male_4>ind[i] & pleio_evo$tau_male_4<ind[i+1]]=i
}

boxplot(log(as.numeric(paste(pleio_evo$F_value)))~group_tau,names=seq(0.15,0.95,0.2)
        ,xlab = "tissue specificity (group median)",ylab = "ln(F-value)")


anc_var_1=data.frame(FBgn=row.names(dat_B_H_all),anc_v=apply(dat_B_WY_old,1,var))
pleio_evo=merge(pleio_evo,anc_var_1,all=T,by="FBgn")
anc_var_2=data.frame(FBgn=row.names(dat_B_H_all),anc_v=apply(dat_B_WY_new,1,var))
pleio_evo=merge(pleio_evo,anc_var_2,all=T,by="FBgn")

fc_tau$avg_fc=apply(fc_tau[,2:31],1,mean)

pleio_evo=merge(pleio_evo,fc_tau[,c(1,34)],by="FBgn")

pleio_evo=pleio_evo[,c(1,2,6,9,10,12)]
pleio_evo_sort=pleio_evo[,c(1,3,4,5,6,2)]
pleio_evo_sort_use=na.omit(pleio_evo_sort)

cor(pleio_evo_sort_use$tau_male_4,pleio_evo_sort_use$conn,method = "spearman")

png(filename = "pleio_hetero.png",height = 15,width = 20,units = "cm",res = 600,pointsize = 9)
par(mfrow=c(1,2))
boxplot(log(as.numeric(paste(pleio_evo$F_value)))~group_tau,names=seq(0.15,0.95,0.2)
        ,xlab = "tissue specificity (group median)",ylab = "ln(F-value)",ylim=c(-2,4))
boxplot(log(as.numeric(paste(pleio_evo$F_value)))~group_conn,ylim=c(-2,4),
        names=seq(1,9,2),xlab = "ln(connectivity) (group median)",
        ylab = "ln(F-value)",ylim=c(-2,4))
dev.off()

cor.test(log(as.numeric(paste(pleio_evo_sort_use$F_value))),pleio_evo_sort_use$tau_male_4,method = "spearman")
cor.test(log(as.numeric(paste(pleio_evo_sort_use$F_value))),pleio_evo_sort_use$conn,method = "spearman")

png(filename = "pleio_anc.png",height = 15,width = 20,units = "cm",res = 600,pointsize = 9)
par(mfrow=c(1,2))
boxplot(as.numeric(paste(pleio_evo$anc_v))~group_tau,names=seq(0.15,0.95,0.2),log="y"
        ,xlab = "tissue specificity (group median)",ylab = "phenotypic variance")
boxplot(as.numeric(paste(pleio_evo$anc_v))~group_conn,
        names=seq(1,9,2),xlab = "ln(connectivity) (group median)",log="y"
        ,ylab = "phenotypic variance")
dev.off()

cor.test(log(as.numeric(paste(pleio_evo_sort_use$anc_v))),pleio_evo_sort_use$tau_male_4,method = "spearman")
cor.test(log(as.numeric(paste(pleio_evo_sort_use$anc_v))),pleio_evo_sort_use$conn,method = "spearman")

png(filename = "pleio_fc.png",height = 15,width = 20,units = "cm",res = 600,pointsize = 9)
par(mfrow=c(1,2))
boxplot(log(abs(as.numeric(paste(pleio_evo$avg_fc))))~group_tau,names=seq(0.15,0.95,0.2)
        ,xlab = "tissue specificity (group median)",ylab = "log2FC")
boxplot(log(abs(as.numeric(paste(pleio_evo$avg_fc))))~group_conn,
        names=seq(1,9,2),xlab = "ln(connectivity) (group median)"
        ,ylab = "log2FC")
dev.off()


cor.test(log(as.numeric(paste(pleio_evo_sort_use$avg_fc))),pleio_evo_sort_use$tau_male_4,method = "spearman")
cor.test(log(as.numeric(paste(pleio_evo_sort_use$avg_fc))),pleio_evo_sort_use$conn,method = "spearman")


cor.test(as.numeric(paste(pleio_evo_sort_use$F_value)),pleio_evo_sort_use$anc_v,method = "spearman")

gray(seq(0,1,0.2),alpha = 0.2)
plot(pleio_evo_sort_use$anc_v,as.numeric(paste(pleio_evo_sort_use$F_value)),log="xy",col=gray(seq(0,1,0.2))[group_conn],pch=1)



png(filename = "anc_hetero.png",height = 15,width = 20,units = "cm",res = 600,pointsize = 9)
plot(log(as.numeric(paste(pleio_evo$F_value))), log(as.numeric(paste(pleio_evo$anc_v.x))))
boxplot(log(as.numeric(paste(pleio_evo$F_value)))~group_tau,names=seq(0.15,0.95,0.2)
        ,xlab = "tissue specificity (group median)",ylab = "ln(F-value)",ylim=c(-2,4))
boxplot(log(as.numeric(paste(pleio_evo$F_value)))~group_conn,ylim=c(-2,4),
        names=seq(1,9,2),xlab = "ln(connectivity) (group median)",
        ylab = "ln(F-value)",ylim=c(-2,4))
dev.off()
pleio_evo_sort_use$tau_male_4=1-pleio_evo_sort_use$tau_male_4

write.table(pleio_evo_sort_use,file = "pleio_evo.txt",sep = "\t",col.names = T,row.names = F, quote = F)

pleio_evo_sort_use = read.table(file = "pleio_evo.txt",sep = "\t",header = T,stringsAsFactors = F)
row.names(pleio_evo_sort_use)=pleio_evo_sort_use[,1]

#repr_gene = genesInTerm(tgd1,"GO:0032504")[[1]]
#pleio_evo_sort_use_no_repr_gene = pleio_evo_sort_use[!row.names(pleio_evo_sort_use)%in%repr_gene,]


cor.test(log(as.numeric(paste(pleio_evo_sort_use_no_repr_gene$F_value))),pleio_evo_sort_use_no_repr_gene$tau_male_4,method = "spearman")
cor.test(log(as.numeric(paste(pleio_evo_sort_use_no_repr_gene$F_value))),pleio_evo_sort_use_no_repr_gene$conn,method = "spearman")

cor.test(log(as.numeric(paste(pleio_evo_sort_use_no_repr_gene$anc_v))),pleio_evo_sort_use_no_repr_gene$tau_male_4,method = "spearman")
cor.test(log(as.numeric(paste(pleio_evo_sort_use_no_repr_gene$anc_v))),pleio_evo_sort_use_no_repr_gene$conn,method = "spearman")

cor.test(log(as.numeric(paste(pleio_evo_sort_use_no_repr_gene$anc_v))),log(as.numeric(paste(pleio_evo_sort_use_no_repr_gene$F_value))),method = "spearman")

plot(as.numeric(paste(pleio_evo_sort_use_no_repr_gene$anc_v)),as.numeric(paste(pleio_evo_sort_use_no_repr_gene$F_value)))




png(filename = "cor_anc_heter.png",height = 15,width = 10,units = "cm",res = 600,pointsize = 9)
boxplot(log(as.numeric(paste(pleio_evo_sort_use$F_value)))~group_anc,ylim=c(-2,4),ylab="Heterogeneity",
        xlab="ancestral variance",names=c("0~20%","20~40%","40~60%","60~80%","80~100%"))
dev.off()


ind=quantile(pleio_evo_sort_use$tau_male_4,seq(0,1,0.1))
group_tau=c()
for (i in 1:10) {
  group_tau[pleio_evo_sort_use$tau_male_4>=ind[i] & pleio_evo_sort_use$tau_male_4<ind[i+1]]=i
}

ind=quantile(pleio_evo_sort_use$conn,seq(0,1,0.1))
group_conn=c()
for (i in 1:10) {
  group_conn[pleio_evo_sort_use$conn>=ind[i] & pleio_evo_sort_use$conn<ind[i+1]]=i
}


dat.use=tapply(log(as.numeric(paste(pleio_evo_sort_use$F_value))), group_anc, function(x) c(mean(x),sd(x)/sqrt(length(x))))

plot(1,1,xlim = c(1,5),ylim = c(0.4,0.9),xlab = "expression variance in ancestral population",
    ylab = "heterogeneity",xaxt= "n")
axis(1,at=1:5,labels = c("20%","40%","60%","80%","100%"))

for (i in 1:5) {
  arrows(x0 = i,x1 = i,y0 = dat.use[[i]][1]-dat.use[[i]][2],y1 = dat.use[[i]][1]+dat.use[[i]][2],
        angle = 90,code = 3,length = 0.1)
  points(i,dat.use[[i]][1],pch=19)
}




cor.test(log(as.numeric(paste(pleio_evo_sort_use$F_value))),pleio_evo_sort_use$tau_male_4,method = "spearman")
cor.test(log(as.numeric(paste(pleio_evo_sort_use$F_value))),pleio_evo_sort_use$conn,method = "spearman")




plot(log2(abs(pleio_evo_sort_use$avg_fc)),log(pleio_evo_sort_use$F_value),
     xlab = "absolute evolutionary change",xlim = c(-10,1),
     ylab = "heterogeneity",xaxt="n")
axis(1,at = c(-8,-4,0),labels = round(2^c(-8,-4,0),3))
text(-10,6,labels = "rho=0.04",pos = 4,col = "red")

cor.test(log(abs(pleio_evo_sort_use$avg_fc)),log(pleio_evo_sort_use$F_value),method = "kendall")





dat_use_filtered=log(dat_use[row.names(pleio_evo_sort_use),])
pleio_evo_sort_use$avg_expr

plot(pleio_evo_sort_use$avg_expr,abs(pleio_evo_sort_use$avg_fc),ylim = c(0,3),
     ylab = "absolute evolutionary response",
     xlab = "avg. expression")
text(2,3,labels = "rho= - 0.05",col = "red")

plot(pleio_evo_sort_use$avg_expr,log(pleio_evo_sort_use$F_value),
     ylab = "heterogeneity",
     xlab = "avg. expression")
text(8,6,labels = "rho = 0.1",col = "red")

cor(pleio_evo_sort_use$avg_expr,log(pleio_evo_sort_use$F_value),method = "spearman")
cor(pleio_evo_sort_use$avg_expr,pleio_evo_sort_use$anc_v.x,method = "pearson")

cor(apply(dat_H_WY_old_v,1,var)[row.names(pleio_evo_sort_use)],log(pleio_evo_sort_use$F_value),method = "spearman")

dat_mean=data.frame(FBgn=names(apply(dat_B_H_filtered, 1, mean)),mean=apply(log(cpm(dat_B_H_filtered)), 1, mean))
pleio_evo_sort_add_mean=merge(dat_mean,pleio_evo_sort_use,by="FBgn")
pleio_evo_sort_add_remove0=pleio_evo_sort_add_mean[-which(pleio_evo_sort_add_mean$conn==0),]

colnames(pleio_evo_sort_add_remove0)=c("FBgn","mean","tau","connectivity","anc_v","fc","f-value")
write.table(pleio_evo_sort_add_remove0,file = "pleio_evo_v2.txt",sep = "\t",col.names = T,row.names = F, quote = F)

pleio_use

