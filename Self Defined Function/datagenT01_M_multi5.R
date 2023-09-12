####Function used to generate recurrent events M_Q1 and M_Q3
M_pre<-function(M_ds){
  M_ds=M_ds[order(M_ds$ID,M_ds$stoptime),]
  M_ds=M_ds[!duplicated( M_ds[ , c("ID","stoptime","event","event_Tp")]),]
}

datagenT01_M_multi<-function(X,X2,T0,T1,vi,tau,
                             lambda.m1,lambda.m2,lambda.m3,lambda.m4,lambda.m5,
                             t.jump.m1,t.jump.m2,t.jump.m3,t.jump.m4,t.jump.m5,
                             betaz1,betaz2,betaz3,betaz4,betaz5,
                             betax1,betax2,betax3,betax4,betax5,
                             delta_X1,delta_X2,delta_X3,delta_X4,delta_X5,event_order,event_z){
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all1=c(betaz1,betax1,delta_X1)
  risk1_0=c(exp(cbind(T0,X,vi)%*%beta.all1))
  risk1_1=c(exp(cbind(T1,X,vi)%*%beta.all1))
  
  beta.all2=c(betaz2,betax2,delta_X2)
  risk2_0=c(exp(cbind(T0,X,vi)%*%beta.all2))
  risk2_1=c(exp(cbind(T1,X,vi)%*%beta.all2))
  
  beta.all3=c(betaz3,betax3,delta_X3)
  risk3_0=c(exp(cbind(T0,X,vi)%*%beta.all3))
  risk3_1=c(exp(cbind(T1,X,vi)%*%beta.all3))
  
  beta.all4=c(betaz4,betax4,delta_X4)
  risk4_0=c(exp(cbind(T0,X2,vi)%*%beta.all4))
  risk4_1=c(exp(cbind(T1,X2,vi)%*%beta.all4))
  
  beta.all5=c(betaz5,betax5,delta_X5)
  risk5_0=c(exp(cbind(T0,X2,vi)%*%beta.all5))
  risk5_1=c(exp(cbind(T1,X2,vi)%*%beta.all5))
  
  ####generate survival time until all individual has follow-up at lest time tau
  M1_0=mysimRec(t.jump.m1,lambda.m1,tau,risk1_0,recurrent=TRUE)
  M1_1=mysimRec(t.jump.m1,lambda.m1,tau,risk1_1,recurrent=TRUE)
  
  M2_0=mysimRec(t.jump.m2,lambda.m2,tau,risk2_0,recurrent=TRUE)
  M2_1=mysimRec(t.jump.m2,lambda.m2,tau,risk2_1,recurrent=TRUE)
  
  M3_0=mysimRec(t.jump.m3,lambda.m3,tau,risk3_0,recurrent=TRUE)
  M3_1=mysimRec(t.jump.m3,lambda.m3,tau,risk3_1,recurrent=TRUE)
  
  M4_0=mysimRec(t.jump.m4,lambda.m4,tau,risk4_0,recurrent=TRUE)
  M4_1=mysimRec(t.jump.m4,lambda.m4,tau,risk4_1,recurrent=TRUE)
  
  M5_0=mysimRec(t.jump.m5,lambda.m5,tau,risk5_0,recurrent=TRUE)
  M5_1=mysimRec(t.jump.m5,lambda.m5,tau,risk5_1,recurrent=TRUE)
  
  M1_0$event_Tp=1
  M1_1$event_Tp=1
  M2_0$event_Tp=2
  M2_1$event_Tp=2
  M3_0$event_Tp=3
  M3_1$event_Tp=3
  M4_0$event_Tp=4
  M4_1$event_Tp=4
  M5_0$event_Tp=5
  M5_1$event_Tp=5
  
  M1_0$M=0
  M1_1$M=1
  M2_0$M=0
  M2_1$M=1
  M3_0$M=0
  M3_1$M=1
  M4_0$M=0
  M4_1$M=1
  M5_0$M=0
  M5_1$M=1
  
  M1_0[which(M1_0$event==0),]$event_Tp=0
  M1_1[which(M1_1$event==0),]$event_Tp=0
  M2_0[which(M2_0$event==0),]$event_Tp=0
  M2_1[which(M2_1$event==0),]$event_Tp=0
  M3_0[which(M3_0$event==0),]$event_Tp=0
  M3_1[which(M3_1$event==0),]$event_Tp=0
  M4_0[which(M4_0$event==0),]$event_Tp=0
  M4_1[which(M4_1$event==0),]$event_Tp=0
  M5_0[which(M5_0$event==0),]$event_Tp=0
  M5_1[which(M5_1$event==0),]$event_Tp=0
  
  All_M0=list(M1_0,M2_0,M3_0,M4_0,M5_0)
  All_M1=list(M1_1,M2_1,M3_1,M4_1,M5_1)
  
 
    if (event_z==T0){M_T=All_M0[[event_order]]}
    if (event_z==T1){M_T=All_M1[[event_order]]}
    M_other0=All_M0[-event_order]
    M_other1=All_M1[-event_order]

  
  M_1=rbind(M_T,M_other0[[1]],M_other0[[2]],M_other0[[3]],M_other0[[4]])
  M_2=rbind(M_T,M_other1[[1]],M_other0[[2]],M_other0[[3]],M_other0[[4]])
  M_3=rbind(M_T,M_other0[[1]],M_other1[[2]],M_other0[[3]],M_other0[[4]])
  M_4=rbind(M_T,M_other0[[1]],M_other0[[2]],M_other1[[3]],M_other0[[4]])
  M_5=rbind(M_T,M_other0[[1]],M_other0[[2]],M_other0[[3]],M_other1[[4]])
  
  M_6=rbind(M_T,M_other1[[1]],M_other1[[2]],M_other0[[3]],M_other0[[4]])
  M_7=rbind(M_T,M_other1[[1]],M_other0[[2]],M_other1[[3]],M_other0[[4]])
  M_8=rbind(M_T,M_other1[[1]],M_other0[[2]],M_other0[[3]],M_other1[[4]])
  M_9=rbind(M_T,M_other0[[1]],M_other1[[2]],M_other1[[3]],M_other0[[4]])
  M_10=rbind(M_T,M_other0[[1]],M_other1[[2]],M_other0[[3]],M_other1[[4]])
  M_11=rbind(M_T,M_other0[[1]],M_other0[[2]],M_other1[[3]],M_other1[[4]])
  
  M_12=rbind(M_T,M_other1[[1]],M_other1[[2]],M_other1[[3]],M_other0[[4]])
  M_13=rbind(M_T,M_other1[[1]],M_other1[[2]],M_other0[[3]],M_other1[[4]])
  M_14=rbind(M_T,M_other1[[1]],M_other0[[2]],M_other1[[3]],M_other1[[4]])
  M_15=rbind(M_T,M_other0[[1]],M_other1[[2]],M_other1[[3]],M_other1[[4]])
  M_16=rbind(M_T,M_other1[[1]],M_other1[[2]],M_other1[[3]],M_other1[[4]])


  
  M_1=M_pre(M_1)
  M_2=M_pre(M_2)
  M_3=M_pre(M_3)
  M_4=M_pre(M_4)
  M_5=M_pre(M_5)
  M_6=M_pre(M_6)
  M_7=M_pre(M_7)
  M_8=M_pre(M_8)
  M_9=M_pre(M_9)
  M_10=M_pre(M_10)
  M_11=M_pre(M_11)
  M_12=M_pre(M_12)
  M_13=M_pre(M_13)
  M_14=M_pre(M_14)
  M_15=M_pre(M_15)
  M_16=M_pre(M_16)
  
  return(list( 
    M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,M_13,M_14,M_15,M_16))
 
}

