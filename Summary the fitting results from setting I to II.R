
####################################################################
##########################True parameters###########################
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
shape=1.5
scale=1

################Computing log() for lambda.m, lambda.y and sigmav###############

log_lam_m1_1=log(lambda.m1[1])
log_lam_m1_2=log(lambda.m1[2])
log_lam_m1_3=log(lambda.m1[3])
log_lam_m1_4=log(lambda.m1[4])

log_lam_m2_1=log(lambda.m2[1])
log_lam_m2_2=log(lambda.m2[2])
log_lam_m2_3=log(lambda.m2[3])
log_lam_m2_4=log(lambda.m2[4])

log_lam_y1=log(lambda.y[1])
log_lam_y2=log(lambda.y[2])
log_lam_y3=log(lambda.y[3])
log_lam_y4=log(lambda.y[4])



##############################################################################
###################################Setting I###################################
True_est=cbind(log_lam_m1_1,log_lam_m1_2,log_lam_m1_3,log_lam_m1_4,
               log_lam_m2_1,log_lam_m2_2,log_lam_m2_3,log_lam_m2_4,
               log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
               betax1,betax2,betaz1,betaz2,
               deltax1,deltax2,etax, etaz,
               delta1, etam1, etam2)

#####Import the estimation and covariance######


est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  est.next=read.table(paste("estImulti",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covImulti",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}
names(est)<-c("log_lam_m1_1","log_lam_m1_2","log_lam_m1_3","log_lam_m1_4",
              "log_lam_m2_1","log_lam_m2_2","log_lam_m2_3","log_lam_m2_4",
              "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
              "betax1","betax2","betaz1","betaz2",
              "deltax1","deltax2","etax", "etaz",
              "delta1", "etam1", "etam2") 
names(var)<-c("log_lam_m1_1","log_lam_m1_2","log_lam_m1_3","log_lam_m1_4",
              "log_lam_m2_1","log_lam_m2_2","log_lam_m2_3","log_lam_m2_4",
              "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
              "betax1","betax2","betaz1","betaz2",
              "deltax1","deltax2","etax", "etaz",
              "delta1", "etam1", "etam2") 


est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim1_estbiaset2",".csv",sep=""),row.names = TRUE)


##############################################################################
###################################Setting II###################################
True_est=cbind(log_lam_m1_1,log_lam_m1_2,log_lam_m1_3,log_lam_m1_4,
               log_lam_m2_1,log_lam_m2_2,log_lam_m2_3,log_lam_m2_4,
               log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
               betax1,betax2,betaz1,betaz2,
               deltax1,deltax2,etax, etaz,
               delta1, etam1, etam2)

#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  est.next=read.table(paste("estIImulti",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covIImulti",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}
names(est)<-c("log_lam_m1_1","log_lam_m1_2","log_lam_m1_3","log_lam_m1_4",
              "log_lam_m2_1","log_lam_m2_2","log_lam_m2_3","log_lam_m2_4",
              "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
              "betax1","betax2","betaz1","betaz2",
              "deltax1","deltax2","etax", "etaz",
              "delta1", "etam1", "etam2") 
names(var)<-c("log_lam_m1_1","log_lam_m1_2","log_lam_m1_3","log_lam_m1_4",
              "log_lam_m2_1","log_lam_m2_2","log_lam_m2_3","log_lam_m2_4",
              "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
              "betax1","betax2","betaz1","betaz2",
              "deltax1","deltax2","etax", "etaz",
              "delta1", "etam1", "etam2") 


est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim2_estbiaset2",".csv",sep=""),row.names = TRUE)


