
library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
###load in self-defined functions#####
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagenT01_M_multi.R")
source("S.R")
source("mymed.R")
source("datagen_X.R")

#######################compute true effects using large data#####################################
#################################################################################################
##### TRUE parameter ######
n=10000
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

est=cbind(log1_r1=log(lambda.m1[1]),log1_r2=log(lambda.m1[2]),log1_r3=log(lambda.m1[3]),log1_r4=log(lambda.m1[4]),
          log2_r1=log(lambda.m2[1]),log2_r2=log(lambda.m2[2]),log2_r3=log(lambda.m2[3]),log2_r4=log(lambda.m2[4]),
          log_h1=log(lambda.y[1]),log_h2=log(lambda.y[2]),log_h3=log(lambda.y[3]),log_h4=log(lambda.y[4]),
          betax1=betax1,betax2=betax2,
          betaz1=betaz1,betaz2=betaz2,
          deltax1=deltax1,deltax2=deltax2,
          etax=etax,etaz=etaz,delta1=delta1,etam1=etam1,etam2=etam2)
est=as.data.frame(est)
#####Generate a large dataset#######


set.seed(8888)
vi=rnorm(n)*sigmav
data=datagen_X(n,vi,tau,betaz1,betaz2,betax1,betax2,lambda.m1,lambda.m2,t.jump.m1,t.jump.m2,
               etaz,etax,etam1,etam2,delta1,deltax1,deltax2,delta2_1,delta2_2,
               lambda.y,t.jump.y,cen=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
tseq=seq(1,8,by=1)


#####Get the true value####

truevalue=mymed(X,n,est,tseq,B=1000)
truevalue_in=c(truevalue$NDE,truevalue$NIE,truevalue$NIE1,truevalue$NIE2)
truevalue1=truevalue_in$x
###handle each simulated dataset

allest=NULL
allsd=NULL
myest=NULL

for (i in 1:200){
  myBest=NULL
  myest1=read.csv(paste("est1_multi",as.character(i),".csv",sep=""))
	myest=c(myest1$NDE,myest1$NIE,myest1$NIE1, myest1$NIE2)
	for (batch in 1:5){
	  skip_to_next <- FALSE
	  tryCatch({myBest_temp=read.csv(paste("Boot1_multilres",as.character(i),"batch=",as.character(batch),".csv",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
	  
	  if(skip_to_next) { next }
	  if (!is.na(sum(myBest_temp))){
	  myBest=rbind(myBest,myBest_temp)
	  }
	  
	}
	
	if (!is.null(myBest)){
	myBest[abs(myBest)>1]<-NA
	myest.sd=apply(myBest,2,sd,na.rm=TRUE)
	allest=rbind(allest,myest)
	allsd=rbind(allsd,myest.sd)
	}
}

bias=apply(allest,2,mean,na.rm=TRUE)-truevalue1
ese=apply(allsd*is.finite(allsd),2,mean,na.rm=TRUE)
meSE=apply(allsd,2,median, na.rm=TRUE)
sd=apply(allest,2,sd,na.rm=TRUE)
cr=apply((abs(allest-rep(1,nrow(allest))%o%truevalue1)/allsd)<1.96,2,mean,na.rm=TRUE)

###might change file name to make it clear the parameter setting
write.csv(rbind(truevalue1,bias,ese,sd,meSE,cr),"sim1_multi_NIENDE_sum.csv")



