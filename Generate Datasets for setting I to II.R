##Load the package

library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
###load in all self defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagen_X.R")


################################################################################################
####################################Setting###################################################
n=200
k=1
tau=10
t.jump.m1=c(2,4,6)
t.jump.m2=c(2,4,6)
t.jump.y=c(2,4,6)
lambda.m1=c(0.65,0.5,0.45,0.35)
lambda.m2=c(0.5,0.35,0.3,0.25)
lambda.y=c(0.033,0.084,0.1,0.6)

betax1=rep(0.2,k)
betax2=rep(0.15,k)
betaz1=0.35
betaz2=0.4
etaz=0.35
etax=rep(0.15,k)
etam1=0.18
etam2=0.15
delta1=0.7
deltax1=0.7
deltax2=0.84
delta2_1=0
delta2_2=0

sigmav=1
shape=1.5
scale=1

#################################Generate 200 datasets from Setting I#########################
delta2_1=0
delta2_2=0

for (iter in 1:200){
  set.seed(8888+iter)
  vi=rnorm(n)*sigmav
  data=datagen_X(n,vi,tau,betaz1,betaz2,betax1,betax2,lambda.m1,lambda.m2,t.jump.m1,t.jump.m2,
                 etaz,etax,etam1,etam2,delta1,deltax1,deltax2,delta2_1,delta2_2,
                 lambda.y,t.jump.y,cen=TRUE)
  
  write.csv(data$alldata,paste("sim1_multi",as.character(iter),".csv",sep=""),row.names=F)
}

#################################Generate 200 datasets from Setting II#########################
for (iter in 1:200){
  set.seed(8888+iter)
  v_gamma=rgamma(n,shape = 1.5, scale = 1)
  vi=log(v_gamma)
  data=datagen_X(n,vi,tau,betaz1,betaz2,betax1,betax2,lambda.m1,lambda.m2,t.jump.m1,t.jump.m2,
                 etaz,etax,etam1,etam2,delta1,deltax1,deltax2,delta2_1,delta2_2,
                 lambda.y,t.jump.y,cen=TRUE)
  write.csv(data$alldata,paste("sim2_multi",as.character(iter),".csv",sep=""),row.names=F)
}

