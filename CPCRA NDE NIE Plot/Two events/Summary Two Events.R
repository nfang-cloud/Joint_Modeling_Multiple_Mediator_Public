
library(sas7bdat)
library(dplyr)
library(survival)

MM=20
tseq=seq(0,15*30,by=30)

#########################################################
##########First summary the real estimates of NIE/NDE####
#########################################################
###################CD4##########################
resCD4_NDE=matrix(data=NA,nrow=MM,ncol=length(tseq))
resCD4_NIE=matrix(data=NA,nrow=MM,ncol=length(tseq))
resCD4_NIE1=matrix(data=NA,nrow=MM,ncol=length(tseq))
resCD4_NIE2=matrix(data=NA,nrow=MM,ncol=length(tseq))
resCD4_TE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_CD4=NULL
for (mm in 1:MM){
  res_CD4=read.csv(paste("Res_Multi_CD4",as.character(1111),"batch=",as.character(mm-1),".csv",sep=""))
  resCD4_NDE[mm,]=c(0,res_CD4$NDE)
  
  resCD4_NIE[mm,]=c(0,res_CD4$NIE)
  resCD4_NIE1[mm,]=c(0,res_CD4$NIE1)
  resCD4_NIE2[mm,]=c(0,res_CD4$NIE2)
  resCD4_TE[mm,]=c(0,res_CD4$TE)
  
}

RCD4_NDE=apply(resCD4_NDE,2,mean,na.rm=TRUE)
RCD4_NIE=apply(resCD4_NIE,2,mean,na.rm=TRUE)
RCD4_NIE1=apply(resCD4_NIE1,2,mean,na.rm=TRUE)
RCD4_NIE2=apply(resCD4_NIE2,2,mean,na.rm=TRUE)
RCD4_TE=apply(resCD4_TE,2,mean,na.rm=TRUE)


################################################
##########Summary the bootstrap of NIE/NDE####
################################################
######average the batch####
Ba_tmp=20
simseed=1000
tmpresCD4_NDE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpresCD4_NIE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpresCD4_NIE1=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpresCD4_NIE2=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpresCD4_TE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
AresCD4_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
AresCD4_NIE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
AresCD4_NIE1=matrix(data=NA,nrow=simseed,ncol=length(tseq))
AresCD4_NIE2=matrix(data=NA,nrow=simseed,ncol=length(tseq))
AresCD4_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresCD4_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresCD4_NIE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresCD4_NIE1=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresCD4_NIE2=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresCD4_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))


BresCD4_tmp=NULL
for (sim in 1:simseed){
  for (mm in 1:Ba_tmp){
    BresCD4_tmp=read.csv(paste("Bres_Multi_CD4",as.character(1111+sim-1),"batch=",as.character(mm),".csv",sep=""))
    tmpresCD4_NDE[mm,]=c(0,BresCD4_tmp$NDE)
    tmpresCD4_NIE[mm,]=c(0,BresCD4_tmp$NIE)
    tmpresCD4_NIE1[mm,]=c(0,BresCD4_tmp$NIE1)
    tmpresCD4_NIE2[mm,]=c(0,BresCD4_tmp$NIE2)
    tmpresCD4_TE[mm,]=c(0,BresCD4_tmp$TE)
    
  }
  AresCD4_NDE[sim,]=apply(tmpresCD4_NDE,2,mean,na.rm=TRUE)
  AresCD4_NIE[sim,]=apply(tmpresCD4_NIE,2,mean,na.rm=TRUE)
  AresCD4_NIE1[sim,]=apply(tmpresCD4_NIE1,2,mean,na.rm=TRUE)
  AresCD4_NIE2[sim,]=apply(tmpresCD4_NIE2,2,mean,na.rm=TRUE)
  AresCD4_TE[sim,]=apply(tmpresCD4_TE,2,mean,na.rm=TRUE)
  SEresCD4_NDE[sim,]=apply(tmpresCD4_NDE,2,sd,na.rm=TRUE)
  SEresCD4_NIE[sim,]=apply(tmpresCD4_NIE,2,sd,na.rm=TRUE)
  SEresCD4_NIE1[sim,]=apply(tmpresCD4_NIE1,2,sd,na.rm=TRUE)
  SEresCD4_NIE2[sim,]=apply(tmpresCD4_NIE2,2,sd,na.rm=TRUE)
  SEresCD4_TE[sim,]=apply(tmpresCD4_TE,2,sd,na.rm=TRUE)
  
}


SDCD4_NDE=apply(AresCD4_NDE,2,sd,na.rm=TRUE)
SDCD4_NIE=apply(AresCD4_NIE,2,sd,na.rm=TRUE)
SDCD4_NIE1=apply(AresCD4_NIE1,2,sd,na.rm=TRUE)
SDCD4_NIE2=apply(AresCD4_NIE2,2,sd,na.rm=TRUE)
SDCD4_TE=apply(AresCD4_TE,2,sd,na.rm=TRUE)

SENDE_CD4=apply(SEresCD4_NDE,2,mean,na.rm=TRUE)
SENIE_CD4=apply(SEresCD4_NIE,2,mean,na.rm=TRUE)
SENIE1_CD4=apply(SEresCD4_NIE1,2,mean,na.rm=TRUE)
SENIE2_CD4=apply(SEresCD4_NIE2,2,mean,na.rm=TRUE)
SETE_CD4=apply(SEresCD4_TE,2,mean,na.rm=TRUE)

QLCD4_NDE=apply(AresCD4_NDE,2,quantile, probs=0.025,na.rm=TRUE)
QUCD4_NDE=apply(AresCD4_NDE,2,quantile, probs=0.975,na.rm=TRUE)

QLCD4_NIE=apply(AresCD4_NIE,2,quantile, probs=0.025,na.rm=TRUE)
QUCD4_NIE=apply(AresCD4_NIE,2,quantile, probs=0.975,na.rm=TRUE)

QLCD4_NIE1=apply(AresCD4_NIE1,2,quantile, probs=0.025,na.rm=TRUE)
QUCD4_NIE1=apply(AresCD4_NIE1,2,quantile, probs=0.975,na.rm=TRUE)

QLCD4_NIE2=apply(AresCD4_NIE2,2,quantile, probs=0.025,na.rm=TRUE)
QUCD4_NIE2=apply(AresCD4_NIE2,2,quantile, probs=0.975,na.rm=TRUE)

QLCD4_TE=apply(AresCD4_TE,2,quantile, probs=0.025,na.rm=TRUE)
QUCD4_TE=apply(AresCD4_TE,2,quantile, probs=0.975,na.rm=TRUE)

##################################################
#####################Plot########################

##################Alternative methods#################
covdata=read.sas7bdat("covar.sas7bdat")
covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$prevoi=covdata$PREVOI
Z_CD4=covdata[,"CD4BL"]/100
Q1=unname(quantile(Z_CD4,na.rm=TRUE, probs=0.25))
Q3=unname(quantile(Z_CD4,na.rm=TRUE, probs=0.75))

surv_CD4=NULL
surv_CD4=read.sas7bdat("surv.sas7bdat")
surv_CD4$event=surv_CD4$event*(1/2)
surv_CD4$cd4bl=Z_CD4



############Cox#############
surv_CD4=cbind(surv_CD4,covdata)
cox_fit_CD4 <- coxph(Surv(fup*30,event)~cd4bl+trt+gender+hemobl+stratum+prevoi, data=surv_CD4)

data0_CD4=data1_CD4=surv_CD4
data0_CD4$cd4bl=Q1
data1_CD4$cd4bl=Q3

fit2CD4_1<-survfit(cox_fit_CD4,newdata=data1_CD4)
fit2CD4_0<-survfit(cox_fit_CD4,newdata=data0_CD4)
CoxCD4_S1<-summary(fit2CD4_1, time=tseq)
CoxCD4_S0<-summary(fit2CD4_0, time=tseq)
fit2_CD4<-apply(CoxCD4_S1[["surv"]]-CoxCD4_S0[["surv"]],1,mean,na.rm=TRUE)

######Treatment#########################

MM=20
tseq=seq(0,15*30,by=30)

#########################################################
##########First summary the real estimates of NIE/NDE####
#########################################################
res_NDE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_NIE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_NIE1=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_NIE2=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_TE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res=NULL
for (mm in 1:MM){
  res=read.csv(paste("Res_Multi_trt",as.character(1111),"batch=",as.character(mm-1),".csv",sep=""))
  res_NDE[mm,]=c(0,res$NDE)
  
  res_NIE[mm,]=c(0,res$NIE)
  res_NIE1[mm,]=c(0,res$NIE1)
  res_NIE2[mm,]=c(0,res$NIE2)
  res_TE[mm,]=c(0,res$TE)
  
}

R_NDE=apply(res_NDE,2,mean,na.rm=TRUE)
R_NIE=apply(res_NIE,2,mean,na.rm=TRUE)
R_NIE1=apply(res_NIE1,2,mean,na.rm=TRUE)
R_NIE2=apply(res_NIE2,2,mean,na.rm=TRUE)
R_TE=apply(res_TE,2,mean,na.rm=TRUE)


################################################
##########Summary the bootstrap of NIE/NDE####
################################################
######average the batch####
Ba_tmp=20
simseed=1000
tmpres_NDE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_NIE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_NIE1=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_NIE2=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_TE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
Ares_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_NIE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_NIE1=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_NIE2=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NIE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NIE1=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NIE2=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))


Bres_tmp=NULL
for (sim in 1:simseed){
  for (mm in 1:Ba_tmp){
    
    skip_to_next <- FALSE
    tryCatch({ Bres_tmp=read.csv(paste("Bres_Multi_trt",as.character(1111+sim-1),"batch=",as.character(mm),".csv",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    if (!is.na(sum(Bres_tmp))){
      
      tmpres_NDE[mm,]=c(0,Bres_tmp$NDE)
      tmpres_NIE[mm,]=c(0,Bres_tmp$NIE)
      tmpres_NIE1[mm,]=c(0,Bres_tmp$NIE1)
      tmpres_NIE2[mm,]=c(0,Bres_tmp$NIE2)
      tmpres_TE[mm,]=c(0,Bres_tmp$TE)
    }
    
    
  }
  Ares_NDE[sim,]=apply(tmpres_NDE,2,mean,na.rm=TRUE)
  Ares_NIE[sim,]=apply(tmpres_NIE,2,mean,na.rm=TRUE)
  Ares_NIE1[sim,]=apply(tmpres_NIE1,2,mean,na.rm=TRUE)
  Ares_NIE2[sim,]=apply(tmpres_NIE2,2,mean,na.rm=TRUE)
  Ares_TE[sim,]=apply(tmpres_TE,2,mean,na.rm=TRUE)
  SEres_NDE[sim,]=apply(tmpres_NDE,2,sd,na.rm=TRUE)
  SEres_NIE[sim,]=apply(tmpres_NIE,2,sd,na.rm=TRUE)
  SEres_NIE1[sim,]=apply(tmpres_NIE1,2,sd,na.rm=TRUE)
  SEres_NIE2[sim,]=apply(tmpres_NIE2,2,sd,na.rm=TRUE)
  SEres_TE[sim,]=apply(tmpres_TE,2,sd,na.rm=TRUE)
  
}


SD_NDE=apply(Ares_NDE,2,sd,na.rm=TRUE)
SD_NIE=apply(Ares_NIE,2,sd,na.rm=TRUE)
SD_NIE1=apply(Ares_NIE1,2,sd,na.rm=TRUE)
SD_NIE2=apply(Ares_NIE2,2,sd,na.rm=TRUE)
SD_TE=apply(Ares_TE,2,sd,na.rm=TRUE)

SENDE=apply(SEres_NDE,2,mean,na.rm=TRUE)
SENIE=apply(SEres_NIE,2,mean,na.rm=TRUE)
SENIE1=apply(SEres_NIE1,2,mean,na.rm=TRUE)
SENIE2=apply(SEres_NIE2,2,mean,na.rm=TRUE)
SETE=apply(SEres_TE,2,mean,na.rm=TRUE)

QL_NDE=apply(Ares_NDE,2,quantile, probs=0.025,na.rm=TRUE)
QU_NDE=apply(Ares_NDE,2,quantile, probs=0.975,na.rm=TRUE)

QL_NIE=apply(Ares_NIE,2,quantile, probs=0.025,na.rm=TRUE)
QU_NIE=apply(Ares_NIE,2,quantile, probs=0.975,na.rm=TRUE)

QL_NIE1=apply(Ares_NIE1,2,quantile, probs=0.025,na.rm=TRUE)
QU_NIE1=apply(Ares_NIE1,2,quantile, probs=0.975,na.rm=TRUE)

QL_NIE2=apply(Ares_NIE2,2,quantile, probs=0.025,na.rm=TRUE)
QU_NIE2=apply(Ares_NIE2,2,quantile, probs=0.975,na.rm=TRUE)

QL_TE=apply(Ares_TE,2,quantile, probs=0.025,na.rm=TRUE)
QU_TE=apply(Ares_TE,2,quantile, probs=0.975,na.rm=TRUE)

##################################################
#####################Plot########################

##################Alternative methods#################
covdata=read.sas7bdat("covar.sas7bdat")
covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$cd4bl=covdata$CD4BL/100
covdata$prevoi=covdata$PREVOI
Z=covdata[,"trt"]
surv=NULL
surv=read.sas7bdat("surv.sas7bdat")
surv$event=surv$event*(1/2)
surv$trt=Z



############Cox#############
surv=cbind(surv,covdata)
cox_fit <- coxph(Surv(fup*30,event)~cd4bl+trt+gender+hemobl+stratum+prevoi, data=surv)
 
data0=data1=surv
data0$trt=0
data1$trt=1

fit2_1<-survfit(cox_fit,newdata=data1)
fit2_0<-survfit(cox_fit,newdata=data0)
Cox_S1<-summary(fit2_1, time=tseq)
Cox_S0<-summary(fit2_0, time=tseq)
fit2<-apply(Cox_S1[["surv"]]-Cox_S0[["surv"]],1,mean,na.rm=TRUE)

###Plot the NDE/NIE/NIE_k/TE for CD4 and TRT ##

pdf("Effects_Multi2.pdf")
par( mfrow= c(2,2) )
plot(c(R_NDE,R_NIE,R_TE)~c(tseq,tseq,tseq),type="n",ylim=c(-0.10,0.15)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)
lines(R_NDE~tseq,col="red",lty=1,lwd=2)
lines(I(QL_NDE)~tseq,col="red",lty=2,lwd=2)
lines(I(QU_NDE)~tseq,col="red",lty=2,lwd=2)

lines(R_NIE~tseq,col="blue",lty=1,lwd=2)
lines(I(QL_NIE)~tseq,col="blue",lty=2,lwd=2)
lines(I(QU_NIE)~tseq,col="blue",lty=2,lwd=2)

lines(R_TE~tseq,col="black",lty=1,lwd=2)
lines(I(QL_TE)~tseq,col="black",lty=2,lwd=2)
lines(I(QU_TE)~tseq,col="black",lty=2,lwd=2)

lines(fit2~tseq,col="green",lty=1,lwd=2)
title(main = "A. Trt (ddC vs ddI)", 
      cex.main = 1.2,   font.main= 2)
legend("bottomleft",c("NDE","NIE","TE","Cox"),col=c("red","blue","black","green"),lty=1,lwd=2,cex=0.8)
abline(h=0)

plot(c(R_NIE1,R_NIE2)~c(tseq,tseq),type="n",ylim=c(-0.03,0.03)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)

lines(R_NIE1~tseq,col="black",lty=1,lwd=2)
lines(I(QL_NIE1)~tseq,col="black",lty=2,lwd=2)
lines(I(QU_NIE1)~tseq,col="black",lty=2,lwd=2)

lines(R_NIE2~tseq,col="purple",lty=1,lwd=2)
lines(I(QL_NIE2)~tseq,col="purple",lty=2,lwd=2)
lines(I(QU_NIE2)~tseq,col="purple",lty=2,lwd=2)
title(main = "B. Subgroup NIEs Trt (ddC vs ddI)", 
      cex.main = 1.2,   font.main= 2)
legend("bottomleft",c("NIE_Lung","NIE_Other"),col=c("black","purple"),lty=1,lwd=2,cex=0.8)
abline(h=0)



plot(c(RCD4_NDE,RCD4_NIE,RCD4_TE)~c(tseq,tseq,tseq),type="n",ylim=c(-0.15,0.25)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)
lines(RCD4_NDE~tseq,col="red",lty=1,lwd=2)
lines(I(QLCD4_NDE)~tseq,col="red",lty=2,lwd=2)
lines(I(QUCD4_NDE)~tseq,col="red",lty=2,lwd=2)

lines(RCD4_NIE~tseq,col="blue",lty=1,lwd=2)
lines(I(QLCD4_NIE)~tseq,col="blue",lty=2,lwd=2)
lines(I(QUCD4_NIE)~tseq,col="blue",lty=2,lwd=2)

lines(RCD4_TE~tseq,col="black",lty=1,lwd=2)
lines(I(QLCD4_TE)~tseq,col="black",lty=2,lwd=2)
lines(I(QUCD4_TE)~tseq,col="black",lty=2,lwd=2)

lines(fit2_CD4~tseq,col="green",lty=1,lwd=2)
title(main = "C. CD4 (Q3 vs Q1)", 
      cex.main = 1.2,   font.main= 2)
legend("bottomleft",c("NDE","NIE","TE","Cox"),col=c("red","blue","black","green"),lty=1,lwd=2,cex=0.8)
abline(h=0)


plot(c(RCD4_NIE1,RCD4_NIE2)~c(tseq,tseq),type="n",ylim=c(-0.02,0.03)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)

lines(RCD4_NIE1~tseq,col="black",lty=1,lwd=2)
lines(I(QLCD4_NIE1)~tseq,col="black",lty=2,lwd=2)
lines(I(QUCD4_NIE1)~tseq,col="black",lty=2,lwd=2)


lines(RCD4_NIE2~tseq,col="purple",lty=1,lwd=2)
lines(I(QLCD4_NIE2)~tseq,col="purple",lty=2,lwd=2)
lines(I(QUCD4_NIE2)~tseq,col="purple",lty=2,lwd=2)
title(main = "D. Subgroup NIEs CD4 (Q3 vs Q1)", 
      cex.main = 1.2,   font.main= 2)
legend("bottomleft",c("NIE_Lung","NIE_Other"),col=c("black","purple"),lty=1,lwd=2,cex=0.8)
abline(h=0)
dev.off()
