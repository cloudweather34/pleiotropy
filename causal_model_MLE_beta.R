# Rscript of likelihood function implementation and power simulation
# author: Wei-Yun Lai
# date created: 2021.06.01
# date last edited: 2022.08.19
# reviewed by Sheng-Kai Hsu 2022.08.21
# General suggestion: move the liklhood function out and use source to call it for the simulation and real data analysis

setwd("~/Dropbox (PopGen)/backup/Wei-Yun/project2/pleiotropy/")

####likelihood function####
# indirect only: pl -> a -> pa
normal.lik1 = function(y,theta) {
  # data and initial parameter specification
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1] # mean of heterogeneity
  sigma_h = theta[2] # sd of heterogeneity
  mu_a=theta[3] # mean of ancestral variation
  sigma_a = theta[4] # sd of ancestral variation
  rho=theta[5] # correlation coefficent between ancestral variation and parallelism (heterogeneity)
  beta_a=theta[6] # beta: effect size of pleiotropy on ancestral variation 
  n = nrow(y) # sample size
  # log liklelihood in parts
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
  logl=logl_1+logl_2+logl_3+logl_4+logl_5 # sum of log-likelihood
  return(-logl)
}
# direct only: pl -> a & pl -> pa
normal.lik2 = function(y, theta) {
  # data and initial parameter specification
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1]
  mu_a=theta[2]
  sigma_h = theta[3]
  sigma_a = theta[4]
  beta_h=theta[5] # effect of pleiotropy on heterogeneity (parallelism)
  beta_a=theta[6] # effect of pleiotropy on ancestral variation
  n = nrow(y)
  # log-likelihood in parts
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
# both direct and indirect: pl -> a -> pa & pl -> a & pl -> pa
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
# direct + only ancestral effect: pl -> pa & a -> pa
normal.lik4 = function(y, theta) {
  y_obs_h=y[,1]
  y_obs_a=y[,2]
  y_ind=y[,3]
  mu_h=theta[1]
  mu_a=theta[2]
  sigma_h = theta[3]
  sigma_a = theta[4]
  rho = theta[5] # correlation between ancestral variation and parallelism (heterogeneity)
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
# Null: no casual relationship
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


####power analysis####
model_I_logl=matrix(NA,50,5)
model_II_logl=matrix(NA,50,5)
model_III_logl=matrix(NA,50,5)
model_IV_logl=matrix(NA,50,5)
model_V_logl=matrix(NA,50,5)


set.seed(100)
for (i in 1:50) {
  print(i)
  # preset initial parameters (matters little as the ML (maximum likelihood) process optimize parameters)
  theta1=c(1,1,2,2,0.5,1)
  theta2=c(1,1,2,2,1,1)
  theta3=c(1,1,2,2,0.5,1,1)
  theta4=c(1,1,2,2,0.5,1)
  theta5=c(1,1,2,2)
  
  # simulated causal relationship (effect size approximated from real data estimate)
  # model I
  yt1=rep(1:10,each=300)
  yt2=-0.13*yt1+rnorm(3000,0,0.97)
  yt3=0.13*yt2+rnorm(3000,0,0.76)
  y1=cbind(yt3,yt2,yt1)
  ## maximum likelihood
  tmp11=optim(theta1,normal.lik1,y=y1,method = "BFGS")$value
  tmp12=optim(theta2,normal.lik2,y=y1,method = "BFGS")$value
  tmp13=optim(theta3,normal.lik3,y=y1,method = "BFGS")$value
  tmp14=optim(theta4,normal.lik4,y=y1,method = "BFGS")$value
  tmp15=optim(theta5,normal.lik5,y=y1,method = "BFGS")$value
  model_I_logl[i,]=c(tmp11,tmp12,tmp13,tmp14,tmp15)
  
  #model II
  yt1=rep(1:10,each=300)
  yt2=-0.13*yt1+rnorm(3000,0,0.97)
  yt3=-0.04*yt1+rnorm(3000,0,0.76)
  y2=cbind(yt3,yt2,yt1)
  tmp21=optim(theta1,normal.lik1,y=y2,method = "BFGS")$value
  tmp22=optim(theta2,normal.lik2,y=y2,method = "BFGS")$value
  tmp23=optim(theta3,normal.lik3,y=y2,method = "BFGS")$value
  tmp24=optim(theta4,normal.lik4,y=y2,method = "BFGS")$value
  tmp25=optim(theta5,normal.lik5,y=y2,method = "BFGS")$value
  model_II_logl[i,]=c(tmp21,tmp22,tmp23,tmp24,tmp25)
  
  # model III
  yt1=rep(1:10,each=300)
  yt2=-0.13*yt1+rnorm(3000,0,0.97)
  yt3=0.1*yt2+-0.03*yt1+rnorm(3000,0,0.76)
  y3=cbind(yt3,yt2,yt1)
  tmp31=optim(theta1,normal.lik1,y=y3,method = "BFGS")$value
  tmp32=optim(theta2,normal.lik2,y=y3,method = "BFGS")$value
  tmp33=optim(theta3,normal.lik3,y=y3,method = "BFGS")$value
  tmp34=optim(theta4,normal.lik4,y=y3,method = "BFGS")$value
  tmp35=optim(theta5,normal.lik5,y=y3,method = "BFGS")$value
  model_III_logl[i,]=c(tmp31,tmp32,tmp33,tmp34,tmp35)
  
  # model IV
  yt1=rep(1:10,each=300)
  yt2=rnorm(3000,1,0.97)
  yt3=0.1*yt2+-0.03*yt1+rnorm(3000,0,0.76)
  y4=cbind(yt3,yt2,yt1)
  tmp41=optim(theta1,normal.lik1,y=y4,method = "BFGS")$value
  tmp42=optim(theta2,normal.lik2,y=y4,method = "BFGS")$value
  tmp43=optim(theta3,normal.lik3,y=y4,method = "BFGS")$value
  tmp44=optim(theta4,normal.lik4,y=y4,method = "BFGS")$value
  tmp45=optim(theta5,normal.lik5,y=y4,method = "BFGS")$value
  model_IV_logl[i,]=c(tmp41,tmp42,tmp43,tmp44,tmp45)
  
  # model V
  yt1=rep(1:10,each=300)
  yt2=rnorm(3000,1,0.97)
  yt3=rnorm(3000,1,0.76)
  y5=cbind(yt3,yt2,yt1)
  tmp51=optim(theta1,normal.lik1,y=y5,method = "BFGS")$value
  tmp52=optim(theta2,normal.lik2,y=y5,method = "BFGS")$value
  tmp53=optim(theta3,normal.lik3,y=y5,method = "BFGS")$value
  tmp54=optim(theta4,normal.lik4,y=y5,method = "BFGS")$value
  tmp55=optim(theta5,normal.lik5,y=y5,method = "BFGS")$value
  model_V_logl[i,]=c(tmp51,tmp52,tmp53,tmp54,tmp55)
}

# BIC calculation: 2*(-logl)+n*log(N)
n=c(6,6,7,6,4) #number of parameter
model_I_BIC=c()
for (i in 1:5) {
  model_I_BIC=cbind(model_I_BIC,2*model_I_logl[,i]+n[i]*log(3000))
}

model_II_BIC=c()
for (i in 1:5) {
  model_II_BIC=cbind(model_II_BIC,2*model_II_logl[,i]+n[i]*log(3000))
}

model_III_BIC=c()
for (i in 1:5) {
  model_III_BIC=cbind(model_III_BIC,2*model_III_logl[,i]+n[i]*log(3000))
}

model_IV_BIC=c()
for (i in 1:5) {
  model_IV_BIC=cbind(model_IV_BIC,2*model_IV_logl[,i]+n[i]*log(3000))
}

model_V_BIC=c()
for (i in 1:5) {
  model_V_BIC=cbind(model_V_BIC,2*model_V_logl[,i]+n[i]*log(3000))
}


res1=table(apply(model_I_BIC,1,which.min))
res2=table(apply(model_II_BIC,1,which.min))
res3=table(apply(model_III_BIC,1,which.min))
res4=table(apply(model_IV_BIC,1,which.min))
res5=table(apply(model_V_BIC,1,which.min))

#### output figure ####
png(filename = "causal_model_power.png",height = 12,width = 12,units = "cm",res = 600,pointsize = 9)
bp.dat=matrix(c(49,0,1,0,0,
                0,50,0,0,0,
                0,0,50,0,0,
                0,0,0,50,0,
                0,0,0,0,50)/50,5,5) # manual summary of the results (res1..5)
x=barplot(bp.dat,horiz = T,names.arg = c("Model I", "Model II","Model III","Model IV","Model V"),
        ylim=c(0,6),
        col = c(alpha("orange",0.7),alpha("olivedrab2",0.7),alpha("grey",0.7),alpha("darksalmon",0.7),alpha("paleturquoise1",0.7)))
text(0.5,3.1,"100%")
text(0.5,1.9,"100%")
text(0.5,0.7,"98%")
text(0.5,4.25,"100%")
text(0.5,5.5,"100%")
dev.off()
