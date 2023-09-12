
library(MASS)
library(survival)
library(dplyr)
###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("datagenT01_M_multi.R")
source("myprod.R")
source("S.R")
source("mymed_T01_2.R")

####check real data###########
####read in data##############

covdata=read.csv("covdata.csv")
covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$cd4bl=covdata$CD4BL/100
covdata$prevoi=covdata$PREVOI
############Z is CD4BL##########
Z=covdata[,"cd4bl"]
T0=unname(quantile(Z,na.rm=TRUE, probs=0.25))
T1=unname(quantile(Z,na.rm=TRUE, probs=0.75))

tseq=seq(1,15,by=1)

X=as.matrix(covdata[,c("prevoi","gender","stratum","hemobl","trt")])
tau=15
n=nrow(X)

############################Import the Estimation######################
 
est=read.table(paste("estI_OL",".txt",sep=""),header=FALSE,row.names=1)
est = as.data.frame(t(est))
 
t.jump.m1=c(1.9,4.7, 6.6, 9.5)
t.jump.m2=c(2.53,5.05,8.55,12.25)
t.jump.y=c(3.63, 7.3, 10.23, 12.27)
delta2_1=0
delta2_2=0


#######################est1 from the real Estimates###############
 
est1=c(est$logL_r1,est$logL_r2,est$logL_r3,est$logL_r4,est$logL_r5,
        est$logO_r1,est$logO_r2,est$logO_r3,est$logO_r4,est$logO_r5,
        est$log_h1,est$log_h2,est$log_h3,est$log_h4,est$log_h5,
        est$betaL_x1,est$betaL_x2,est$betaL_x3,est$betaL_x4,est$betaL_x5,est$betaL_z,
        est$betaO_x1,est$betaO_x2,est$betaO_x3,est$betaO_x4,est$betaO_x5,est$betaO_z,
        est$deltaL,est$deltaO,
        est$etax1,est$etax2,est$etax3,est$etax4,est$etax5,est$etaz,
        est$etamL,est$etamO,est$delta1)
 
 
 

 
 ##Mediation Analysis##
 B=100000/20
 for (batch in 1:20){
 set.seed(1111)
 Bres1_new=NULL
 Best=est1
 Bsigmav=1
 Blambda.m1=exp(Best[1:5])
 Blambda.m2=exp(Best[6:10])
 Blambda.y=exp(Best[11:15])
 
 Bbetaz1=Best[20]
 Bbetaz2=Best[26]
 Bbetax1=Best[c(16:19,21)]
 Bbetax2=Best[c(22:25,27)]
 Bdelta_X1=Best[28]
 Bdelta_X2=Best[29]
 Betaz=Best[34]
 Betax=Best[c(30:33,35)]
 Betam1=Best[36]
 Betam2=Best[37]
 Bdelta1=Best[38]
 Bdelta2_1=0
 Bdelta2_2=0

 Bres1_new=mymed_T01(X,n,Z,T0,T1,Bsigmav,tau,
                     Bbetaz1,Bbetaz2,Bbetax1,Bbetax2,Blambda.m1,Blambda.m2,
                     t.jump.m1,t.jump.m2,Bdelta_X1,Bdelta_X2,
                     Betaz,Betax,Betam1,Betam2,Bdelta1,Bdelta2_1,Bdelta2_2,
                     Blambda.y,t.jump.y,tseq,B,seed=1111+batch*B)
 
 
 write.csv(Bres1_new,paste("Res_Multi_CD4",as.character(1111),"batch=",as.character(batch),".csv",sep=""),row.names=FALSE)
}
 