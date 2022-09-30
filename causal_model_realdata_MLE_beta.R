# Rscript of causal analysis on real data
# author: Wei-Yun Lai
# date created: 2021.06.22
# date last edited: 2022.08.19
# reviewed by Sheng-Kai Hsu 2022.08.21
# suggestion: 
# 1. use source to call likelihood functions
# 2. is it need to run the causal analysis on scaled data
  #response from WY: did it

rm(list=ls())
setwd("~/Dropbox (PopGen)/backup/Wei-Yun/project2/pleiotropy/")


####likelihood function####
normal.lik1 = function(y,theta) {
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1]
  sigma_h = theta[2]
  mu_a=theta[3]
  sigma_a = theta[4]
  rho=theta[5]
  beta_a=theta[6]
  n = nrow(y)
  logl_1 = 0
  for (i in unique(y_ind)) {
    tmp = sum(y_ind==i)^2/n
    logl_1 = tmp+logl_1
  }
  logl_2 = -.5*n*log(2*pi) -.5*n*log(sigma_h*sigma_h) -.5*n*log(1-rho*rho)
  logl_3 = -(1/(2*sigma_h*sigma_h*(1-rho*rho)))*sum((y_obs_h-mu_h-rho*1/sigma_a*sigma_h*(y_obs_a-mu_a))^2)
  logl_4 =-.5*n*log(2*pi) -.5*n*log(sigma_a*sigma_a) 
  logl_5 = 0
  for (i in unique(y_ind)) {
    tmp= - (1/(2*sigma_a^2))*sum((y_obs_a[which(y_ind==i)]-(mu_a+beta_a*(i-(min(y_ind)+max(y_ind))/2)))**2)
    logl_5=tmp+logl_5  
  }
  logl=logl_1+logl_2+logl_3+logl_4+logl_5
  return(-logl)
}
normal.lik2 = function(y, theta) {
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1]
  mu_a=theta[2]
  sigma_h = theta[3]
  sigma_a = theta[4]
  beta_h=theta[5]
  beta_a=theta[6]
  n = nrow(y)
  logl_1 = 0
  for (i in unique(y_ind)) {
    tmp = sum(y_ind==i)^2/n
    logl_1 = tmp+logl_1
  }
  logl_2 = -.5*n*log(2*pi) -.5*n*log(sigma_a*sigma_a)
  logl_3 = -.5*n*log(2*pi) -.5*n*log(sigma_h*sigma_h)
  logl_4 = 0
  for (i in unique(y_ind)) {
    tmp= - (1/(2*sigma_a^2))*sum((y_obs_a[which(y_ind==i)]-(mu_a+beta_a*(i-(min(y_ind)+max(y_ind))/2)))**2)
    logl_4=tmp+logl_4  
  }
  logl_5 = 0
  for (i in unique(y_ind)) {
    tmp= - (1/(2*sigma_h^2))*sum((y_obs_h[which(y_ind==i)]-(mu_h+beta_h*(i-(min(y_ind)+max(y_ind))/2)))**2)
    logl_5=tmp+logl_5  
  }
  logl=logl_1+logl_2+logl_3+logl_4+logl_5
  return(-logl)
}
normal.lik3 = function(y, theta) {
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1]
  mu_a=theta[2]
  sigma_h = theta[3]
  sigma_a = theta[4]
  rho = theta[5]
  beta_h=theta[6]
  beta_a=theta[7]
  n = nrow(y)
  logl_1 = 0
  for (i in unique(y_ind)) {
    tmp = sum(y_ind==i)^2/n
    logl_1 = tmp+logl_1
  }
  logl_2 = -.5*n*log(2*pi) -.5*n*log(sigma_a*sigma_a)
  logl_3 = -.5*n*log(2*pi) -.5*n*log(sigma_h*sigma_h) -.5*n*log(1-rho*rho)
  logl_4 = 0
  for (i in unique(y_ind)) {
    tmp= - (1/(2*sigma_a^2))*sum((y_obs_a[which(y_ind==i)]-(mu_a+beta_a*(i-(min(y_ind)+max(y_ind))/2)))**2)
    logl_4=tmp+logl_4  
  }
  logl_5 = 0
  for (i in unique(y_ind)) {
    tmp= - 1/(2*sigma_h^2*(1-rho*rho))*sum((y_obs_h[which(y_ind==i)]-(mu_h+beta_h*(i-(min(y_ind)+max(y_ind))/2))-rho*sigma_h/sigma_a*(y_obs_a[which(y_ind==i)]-mu_a))**2)
    logl_5=tmp+logl_5  
  }
  
  logl=logl_1+logl_2+logl_3+logl_4+logl_5
  return(-logl)
}
normal.lik4 = function(y, theta) {
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1]
  mu_a=theta[2]
  sigma_h = theta[3]
  sigma_a = theta[4]
  rho = theta[5]
  beta_h=theta[6]
  n = nrow(y)
  
  logl_1 = 0
  for (i in unique(y_ind)) {
    tmp = sum(y_ind==i)^2/n
    logl_1 = tmp+logl_1
  }
  logl_2 = -.5*n*log(2*pi) -.5*n*log(sigma_a*sigma_a)
  logl_3 = -.5*n*log(2*pi) -.5*n*log(sigma_h*sigma_h) -.5*n*log(1-rho*rho)
  logl_4 = -(1/(2*sigma_a^2))*sum((y_obs_a-mu_a)**2)
  logl_5 = 0
  for (i in unique(y_ind)) {
    tmp= - 1/(2*sigma_h^2*(1-rho*rho))*sum((y_obs_h[which(y_ind==i)]-(mu_h+beta_h*(i-(min(y_ind)+max(y_ind))/2))-rho*sigma_h/sigma_a*(y_obs_a[which(y_ind==i)]-mu_a))**2)
    logl_5=tmp+logl_5  
  }
  
  logl=logl_1+logl_2+logl_3+logl_4+logl_5
  return(-logl)
}
normal.lik5 = function(y, theta) {
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1]
  mu_a=theta[2]
  sigma_h = theta[3]
  sigma_a = theta[4]
  n = nrow(y)
  logl_1 = 0
  for (i in unique(y_ind)) {
    tmp = sum(y_ind==i)^2/n
    logl_1 = tmp+logl_1
  }
  logl_2 = -.5*n*log(2*pi) -.5*n*log(sigma_a*sigma_a)
  logl_3 = -.5*n*log(2*pi) -.5*n*log(sigma_h*sigma_h)
  logl_4 = -(1/(2*sigma_a^2))*sum((y_obs_a-mu_a)**2)
  logl_5 = -(1/(2*sigma_h^2))*sum((y_obs_h-mu_h)**2)
  
  logl=logl_1+logl_2+logl_3+logl_4+logl_5
  return(-logl)
}
####data input####
dat.use=read.table("./pleio_putatively_adaptive_gene.txt",sep = "\t",header = T,stringsAsFactors = F)

#ind=quantile(dat.use$connectivity,seq(0,1,0.1))
#ind[11]=ind[11]+1
ind=c(1,2,4,8,16,32,64,128,256,512,5000)
group_conn=c()
for (i in 1:10) {
  group_conn[dat.use$connectivity>=ind[i] & dat.use$connectivity<ind[i+1]]=i
}

#ind=quantile(dat.use$ts_1.tau,seq(0,1,0.1))
#ind[11]=ind[11]+1
ind=seq(0,1.1,0.1)
group_tau=c()
for (i in 1:10) {
  group_tau[dat.use$ts_1.tau>=ind[i] & dat.use$ts_1.tau<ind[i+1]]=i
}

dat.use$group_conn=scale(group_conn)
dat.use$group_tau=scale(group_tau)
dat.use$anc_v_log=scale(log(dat.use$anc_v))
dat.use$f.value_log=scale(-log(dat.use$f.value))


####likelihood calculation-network connectivity####
n=c(6,6,7,6,4)

theta1=c(1,1,2,2,0.5,1)
theta2=c(1,1,2,2,1,1)
theta3=c(1,1,2,2,0.5,1,1)
theta4=c(1,1,2,2,0.5,1)
theta5=c(1,1,2,2)

model_I_logl=optim(theta1,normal.lik1,y=dat.use[,c(12,11,9)],method = "BFGS")$value
model_II_logl=optim(theta2,normal.lik2,y=dat.use[,c(12,11,9)],method = "BFGS")$value
model_III_logl=optim(theta3,normal.lik3,y=dat.use[,c(12,11,9)],method = "BFGS")$value
model_IV_logl=optim(theta4,normal.lik4,y=dat.use[,c(12,11,9)],method = "BFGS")$value
model_V_logl=optim(theta5,normal.lik5,y=dat.use[,c(12,11,9)],method = "BFGS")$value


#BIC
model_I_BIC=2*model_I_logl+n[1]*log(nrow(dat.use))
model_II_BIC=2*model_II_logl+n[2]*log(nrow(dat.use))
model_III_BIC=2*model_III_logl+n[3]*log(nrow(dat.use))
model_IV_BIC=2*model_IV_logl+n[4]*log(nrow(dat.use))
model_V_BIC=2*model_V_logl+n[5]*log(nrow(dat.use))

print(c(model_I_BIC,model_II_BIC,model_III_BIC,model_IV_BIC,model_V_BIC))


####likelihood calculation-tau####
n=c(6,6,7,6,4)
theta1=c(1,1,2,2,0.5,1)
theta2=c(1,1,2,2,1,1)
theta3=c(1,1,2,2,0.5,1,1)
theta4=c(1,1,2,2,0.5,1)
theta5=c(1,1,2,2)


model_I_logl=optim(theta1,normal.lik1,y=dat.use[,c(12,11,10)],method = "BFGS")$value
model_II_logl=optim(theta2,normal.lik2,y=dat.use[,c(12,11,10)],method = "BFGS")$value
model_III_logl=optim(theta3,normal.lik3,y=dat.use[,c(12,11,10)],method = "BFGS")$value
model_IV_logl=optim(theta4,normal.lik4,y=dat.use[,c(12,11,10)],method = "BFGS")$value
model_V_logl=optim(theta5,normal.lik5,y=dat.use[,c(12,11,10)],method = "BFGS")$value


#BIC
model_I_BIC=2*model_I_logl+n[1]*log(nrow(dat.use))
model_II_BIC=2*model_II_logl+n[2]*log(nrow(dat.use))
model_III_BIC=2*model_III_logl+n[3]*log(nrow(dat.use))
model_IV_BIC=2*model_IV_logl+n[4]*log(nrow(dat.use))
model_V_BIC=2*model_V_logl+n[5]*log(nrow(dat.use))

print(c(model_I_BIC,model_II_BIC,model_III_BIC,model_IV_BIC,model_V_BIC))


####effect size estimation- direct and indirect effect####
# get coefficents from two regression models
a_conn=lm(dat.use$f.value_log~dat.use$group_conn+dat.use$anc_v_log)
b_conn=lm(dat.use$anc_v_log~dat.use$group_conn)
anova(a_conn)

# calculate direct/indirect effect size
direct_conn=a_conn$coefficients[2]
indirect_conn=a_conn$coefficients[3]*b_conn$coefficients[2]
print(c(direct_conn,indirect_conn))

#process repeated for tau
a_tau=lm(dat.use$f.value_log~dat.use$group_tau+dat.use$anc_v_log)
b_tau=lm(dat.use$anc_v_log~dat.use$group_tau)
anova(a_tau)

direct_tau=a_tau$coefficients[2]
indirect_tau=a_tau$coefficients[3]*b_tau$coefficients[2]
print(c(direct_tau,indirect_tau))

####ancestral variation on parallelism####
anc_par=lm(dat.use$f.value_log~dat.use$anc_v_log)
print(anc_par$coefficients[2])

####difference in pleiotropy####
a=lm(dat.use$f.value_log~dat.use$anc_v_log+dat.use$group_tau+dat.use$group_conn)
print(a)
anova(a)
b=lm(dat.use$anc_v_log~dat.use$group_tau+dat.use$group_conn)
print(b)
anova(b)








####boxplot pleiotropy and replicate frequency specturm####

#tissue specificity

png("./supplementary_figure.png",width = 14,height = 7,units = "cm",pointsize = 6,res = 600)
par(mfrow=c(1,2))
low=dat.use$ts_1.tau[dat.use$rep_spe==c(2,3,4)]
inter=dat.use$ts_1.tau[dat.use$rep_spe==c(5,6,7)]
high=dat.use$ts_1.tau[dat.use$rep_spe==c(8,9,10)]
boxplot(low,inter,high,names = c("2~4","5~7","8~10"),ylab="1-tau")

low=log(dat.use$connectivity[dat.use$rep_spe==c(3,4)])
inter=log(dat.use$connectivity[dat.use$rep_spe==c(5,6,7)])
high=log(dat.use$connectivity[dat.use$rep_spe==c(8,9,10)])
boxplot(low,inter,high,names = c("2~4","5~7","8~10"),ylab="network connectivity")
dev.off()


