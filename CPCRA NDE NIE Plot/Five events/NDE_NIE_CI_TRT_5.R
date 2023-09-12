
library(MASS)
library(survival)
library(dplyr)
###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("datagenT01_M_multi5.R")
source("myprod.R")
source("S.R")
source("mymed_T01_5.R")


####check real data###########
####read in data##############
covdata=read.csv("covdata.csv")

covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$cd4bl=covdata$CD4BL/100
covdata$prevoi=covdata$PREVOI
Z=covdata[,"trt"]
T0=0
T1=1

tseq=seq(1,15,by=1)

X=as.matrix(covdata[,c("prevoi","gender","stratum","hemobl","cd4bl")])
X2=as.matrix(covdata[,c("cd4bl")])
tau=15
n=nrow(X)

###Load the estimates from real data analysis####
cov1=read.table("covI_Five.txt",row.names=1) 
est=read.table(paste("estI_Five",".txt",sep=""),header=FALSE,row.names=1)
 est = as.data.frame(t(est))
 
 t.jump.m1=5.5
 t.jump.m2=5.2
 t.jump.m3=6.6
 t.jump.m4=11.65
 t.jump.m5=5.3
 t.jump.y=c(3.63, 7.3, 10.23, 12.27)

 

#######################est1 from the real Estimates###############
 est1=c(est$log_r_PCP1,est$log_r_PCP2,
        est$log_r_MAC1,est$log_r_MAC2,
        est$log_r_CMV1,est$log_r_CMV2,
        est$log_r_WAST1,est$log_r_WAST2,
        est$log_r_TOXO1,est$log_r_TOXO2,
        est$log_h1,est$log_h2,est$log_h3,est$log_h4,est$log_h5,
        est$betaPCP_x1,est$betaPCP_x2,est$betaPCP_x3,est$betaPCP_x4,est$betaPCP_x5,est$betaPCP_z,
        est$betaMAC_x1,est$betaMAC_x2,est$betaMAC_x3,est$betaMAC_x4,est$betaMAC_x5,est$betaMAC_z,
        est$betaCMV_x1,est$betaCMV_x2,est$betaCMV_x3,est$betaCMV_x4,est$betaCMV_x5,est$betaCMV_z,
        est$betaWAST_x5,est$betaWAST_z,
        est$betaTOXO_x5,est$betaTOXO_z,
        est$deltaPCP,est$deltaMAC,est$deltaCMV,est$deltaWAST,est$deltaTOXO,
        est$etax1,est$etax2,est$etax3,est$etax4,est$etax5,est$etaz,
        est$etamPCP,est$etamMAC,est$etamCMV,
        est$etamWAST,est$etamTOXO,est$delta1)
 
 
 ###parametric Bootstrap
 
 B=10000/20
 
 ##Mediation Analysis##
 for (batch in 1:20){
         
        for (simseed in 0:999){
 set.seed(1111+simseed)
 Bres1_new=NULL
 Best=mvrnorm(n = 1, mu=est1, Sigma=cov1)
 Bsigmav=sqrt(1)
 Blambda.m1=exp(Best[1:2])
 Blambda.m2=exp(Best[3:4])
 Blambda.m3=exp(Best[5:6])
 Blambda.m4=exp(Best[7:8])
 Blambda.m5=exp(Best[9:10])
 Blambda.y=exp(Best[11:15])
 
 Bbetaz1=Best[21]
 Bbetaz2=Best[27]
 Bbetaz3=Best[33]
 Bbetaz4=Best[35]
 Bbetaz5=Best[37]
 Bbetax1=Best[c(16:20)]
 Bbetax2=Best[c(22:26)]
 Bbetax3=Best[c(28:32)]
 Bbetax4=Best[34]
 Bbetax5=Best[36]
 Bdelta_X1=Best[38]
 Bdelta_X2=Best[39]
 Bdelta_X3=Best[40]
 Bdelta_X4=Best[41]
 Bdelta_X5=Best[42]
 Betaz=Best[48]
 Betax=Best[c(43:47)]
 Betam1=Best[49]
 Betam2=Best[50]
 Betam3=Best[51]
 Betam4=Best[52]
 Betam5=Best[53]
 Bdelta1=Best[54]

 
 
Bres1_new=mymed_T01(X,X2,n,Z,T0,T1,Bsigmav,tau,
                     Blambda.m1,Blambda.m2,Blambda.m3,Blambda.m4,Blambda.m5,
                     t.jump.m1,t.jump.m2,t.jump.m3,t.jump.m4,t.jump.m5,
                     Blambda.y,t.jump.y,
                     Bbetaz1,Bbetaz2,Bbetaz3,Bbetaz4,Bbetaz5,
                     Bbetax1,Bbetax2,Bbetax3,Bbetax4,Bbetax5,
                     Bdelta_X1,Bdelta_X2,Bdelta_X3,Bdelta_X4,Bdelta_X5,
                     Betaz,Betax,
                     Betam1,Betam2,Betam3,Betam4,Betam5,Bdelta1,
                     tseq,B,seed=1111+batch*B)

 write.csv(Bres1_new[[1]],paste("Bres_MultiFive_trt",as.character(1111+simseed),"batch=",as.character(batch),".csv",sep=""),row.names=FALSE)
 }
}
 
 