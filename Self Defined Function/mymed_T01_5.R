####Mediation analysis for five events####
###Function to compute the product of five events' jump points and lambda###
myprod5<-function(t.jump.risk_1,lambda.risk_1,t.jump.risk_2,lambda.risk_2,
                  t.jump.risk_3,lambda.risk_3,t.jump.risk_4,lambda.risk_4,
                  t.jump.risk_5,lambda.risk_5){
  
      lambda.risk12<-myprod(t.jump.risk_1,lambda.risk_1,t.jump.risk_2,lambda.risk_2)[[2]]
      t.jump.risk12<-myprod(t.jump.risk_1,lambda.risk_1,t.jump.risk_2,lambda.risk_2)[[1]]
      lambda.risk123<-myprod(t.jump.risk12, lambda.risk12,t.jump.risk_3,lambda.risk_3)[[2]]
      t.jump.risk123<-myprod(t.jump.risk12, lambda.risk12,t.jump.risk_3,lambda.risk_3)[[1]]
      lambda.risk1234<-myprod(t.jump.risk123, lambda.risk123,t.jump.risk_4,lambda.risk_4)[[2]]
      t.jump.risk1234<-myprod(t.jump.risk123, lambda.risk123,t.jump.risk_4,lambda.risk_4)[[1]]
      lambda.risk12345<-myprod(t.jump.risk1234, lambda.risk1234,t.jump.risk_5,lambda.risk_5)[[2]]
      t.jump.risk12345<-myprod(t.jump.risk1234, lambda.risk1234,t.jump.risk_5,lambda.risk_5)[[1]]
      return(list(t.jump.risk12345,lambda.risk12345))
}


#####function tom compute the survival probability for subject i
S_prepare<-function(R,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T01,etaz,etax,delta1,tseq){
  Ri=R[which(R$ID==i),]
  R_1i=Ri[which(Ri$event_Tp==1),]
  R_2i=Ri[which(Ri$event_Tp==2),]
  R_3i=Ri[which(Ri$event_Tp==3),]
  R_4i=Ri[which(Ri$event_Tp==4),]
  R_5i=Ri[which(Ri$event_Tp==5),]
  
  lambda.risk_1=exp((0:length(R_1i$stoptime))*(etam1[1]+etam1[2]*vii))
  t.jump.risk_1=R_1i$stoptime
  
  lambda.risk_2=exp((0:length(R_2i$stoptime))*(etam2[1]+etam2[2]*vii))
  t.jump.risk_2=R_2i$stoptime
  
  lambda.risk_3=exp((0:length(R_3i$stoptime))*(etam3[1]+etam3[2]*vii))
  t.jump.risk_3=R_3i$stoptime
  
  lambda.risk_4=exp((0:length(R_4i$stoptime))*(etam4[1]+etam4[2]*vii))
  t.jump.risk_4=R_4i$stoptime
  
  lambda.risk_5=exp((0:length(R_5i$stoptime))*(etam5[1]+etam5[2]*vii))
  t.jump.risk_5=R_5i$stoptime
  
  
  t.jump.risk<-myprod5(t.jump.risk_1,lambda.risk_1,t.jump.risk_2,lambda.risk_2,
                       t.jump.risk_3,lambda.risk_3,t.jump.risk_4,lambda.risk_4,
                       t.jump.risk_5,lambda.risk_5)[[1]]
  lambda.risk<-myprod5(t.jump.risk_1,lambda.risk_1,t.jump.risk_2,lambda.risk_2,
                       t.jump.risk_3,lambda.risk_3,t.jump.risk_4,lambda.risk_4,
                       t.jump.risk_5,lambda.risk_5)[[2]]
  my=myprod(t.jump.y,lambda.y,t.jump.risk,lambda.risk)
  S(Xi,T01,vii,my,etaz,etax,delta1,tseq)
}


####Function to compute the 32 combination of S for preparing sub NIE

NIE_prep<-function(X,X2,n,T0,T1,sigmav,tau,
                  lambda.m1,lambda.m2,lambda.m3,lambda.m4,lambda.m5,
                  t.jump.m1,t.jump.m2,t.jump.m3,t.jump.m4,t.jump.m5,
                  lambda.y,t.jump.y,
                  betaz1,betaz2,betaz3,betaz4,betaz5,
                  betax1,betax2,betax3,betax4,betax5,
                  delta_X1,delta_X2,delta_X3,delta_X4,delta_X5,
                  etaz,etax,
                  etam1,etam2,etam3,etam4,etam5,delta1,
                  event_order,event_z,tseq,B,seed){
  
      avgres=matrix(data=0,nrow=n,ncol=32*length(tseq))  
      avgressq= matrix(data=0,nrow=n,ncol=32*length(tseq))
  
      for (b in 1:B){
        
        set.seed(seed+b)
        vi=rnorm(n,0,sigmav)
        #M_T generates 16 combination of M when event is event_order with treatment event_z
        M_T=datagenT01_M_multi(X,X2,T0,T1,vi,tau,
                              lambda.m1,lambda.m2,lambda.m3,lambda.m4,lambda.m5,
                              t.jump.m1,t.jump.m2,t.jump.m3,t.jump.m4,t.jump.m5,
                              betaz1,betaz2,betaz3,betaz4,betaz5,
                              betax1,betax2,betax3,betax4,betax5,
                              delta_X1,delta_X2,delta_X3,delta_X4,delta_X5,event_order,event_z)
        
        
        M_1=M_T[[1]]
        M_2=M_T[[2]]
        M_3=M_T[[3]]
        M_4=M_T[[4]]
        M_5=M_T[[5]]
        M_6=M_T[[6]]
        M_7=M_T[[7]]
        M_8=M_T[[8]]
        M_9=M_T[[9]]
        M_10=M_T[[10]]
        M_11=M_T[[11]]
        M_12=M_T[[12]]
        M_13=M_T[[13]]
        M_14=M_T[[14]]
        M_15=M_T[[15]]
        M_16=M_T[[16]]
        
        R1=M_1[which(M_1$event==1),][,-c(3,5)]
        R2=M_2[which(M_2$event==1),][,-c(3,5)]
        R3=M_3[which(M_3$event==1),][,-c(3,5)]
        R4=M_4[which(M_4$event==1),][,-c(3,5)]
        R5=M_5[which(M_5$event==1),][,-c(3,5)]
        R6=M_6[which(M_6$event==1),][,-c(3,5)]
        R7=M_7[which(M_7$event==1),][,-c(3,5)]
        R8=M_8[which(M_8$event==1),][,-c(3,5)]
        R9=M_9[which(M_9$event==1),][,-c(3,5)]
        R10=M_10[which(M_10$event==1),][,-c(3,5)]
        R11=M_11[which(M_11$event==1),][,-c(3,5)]
        R12=M_12[which(M_12$event==1),][,-c(3,5)]
        R13=M_13[which(M_13$event==1),][,-c(3,5)]
        R14=M_14[which(M_14$event==1),][,-c(3,5)]
        R15=M_15[which(M_15$event==1),][,-c(3,5)]
        R16=M_16[which(M_16$event==1),][,-c(3,5)]
        
        tmpres=avgres
        etam1=c(etam1,0)
        etam2=c(etam2,0)
        etam3=c(etam3,0)
        etam4=c(etam4,0)
        etam5=c(etam5,0)
        
        #S0_/S1_ here 0/1 is the treatment for survival
        S0_1b=S0_2b=S0_3b=S0_4b=S0_5b=S0_6b=S0_7b=S0_8b=S0_9b=S0_10b=S0_11b=S0_12b=S0_13b=S0_14b=S0_15b=S0_16b=matrix(ncol=length(tseq),nrow=n)
        S1_1b=S1_2b=S1_3b=S1_4b=S1_5b=S1_6b=S1_7b=S1_8b=S1_9b=S1_10b=S1_11b=S1_12b=S1_13b=S1_14b=S1_15b=S1_16b=matrix(ncol=length(tseq),nrow=n)
        
        
        for (i in 1:n){
          Xi=X[i,]
          vii=vi[i]
          
          S0_1b=S_prepare(R1,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_2b=S_prepare(R2,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_3b=S_prepare(R3,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_4b=S_prepare(R4,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_5b=S_prepare(R5,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_6b=S_prepare(R6,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_7b=S_prepare(R7,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_8b=S_prepare(R8,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_9b=S_prepare(R9,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_10b=S_prepare(R10,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_11b=S_prepare(R11,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_12b=S_prepare(R12,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_13b=S_prepare(R13,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_14b=S_prepare(R14,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_15b=S_prepare(R15,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S0_16b=S_prepare(R16,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T0,etaz,etax,delta1,tseq)
          S1_1b=S_prepare(R1,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_2b=S_prepare(R2,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_3b=S_prepare(R3,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_4b=S_prepare(R4,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_5b=S_prepare(R5,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_6b=S_prepare(R6,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_7b=S_prepare(R7,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_8b=S_prepare(R8,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_9b=S_prepare(R9,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_10b=S_prepare(R10,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_11b=S_prepare(R11,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_12b=S_prepare(R12,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_13b=S_prepare(R13,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_14b=S_prepare(R14,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_15b=S_prepare(R15,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          S1_16b=S_prepare(R16,i,etam1,etam2,etam3,etam4,etam5,vii,t.jump.y,lambda.y,Xi,T1,etaz,etax,delta1,tseq)
          
          
          tmpres[i,]= c(S0_1b,S0_2b,S0_3b,S0_4b,S0_5b,S0_6b,S0_7b,S0_8b,
                        S0_9b,S0_10b,S0_11b,S0_12b,S0_13b,S0_14b,S0_15b,S0_16b,
                        S1_1b,S1_2b,S1_3b,S1_4b,S1_5b,S1_6b,S1_7b,S1_8b,
                        S1_9b,S1_10b,S1_11b,S1_12b,S1_13b,S1_14b,S1_15b,S1_16b)
        }
        
        tmpressq=tmpres^2
        avgres=avgres+(tmpres-avgres)/b
        avgressq=avgressq+(tmpressq-avgressq)/b
        
      }
     avgsd=sqrt(avgressq-avgres^2)
     avg_S=apply(avgres,2,mean,na.rm=TRUE)
     sd_S=apply(avgres,2,sd,na.rm=TRUE)

list(avgsd,avg_S,sd_S)
}




mymed_T01<-function(X,X2,n,Z,T0,T1,sigmav,tau,
                    lambda.m1,lambda.m2,lambda.m3,lambda.m4,lambda.m5,
                    t.jump.m1,t.jump.m2,t.jump.m3,t.jump.m4,t.jump.m5,
                    lambda.y,t.jump.y,
                    betaz1,betaz2,betaz3,betaz4,betaz5,
                    betax1,betax2,betax3,betax4,betax5,
                    delta_X1,delta_X2,delta_X3,delta_X4,delta_X5,
                    etaz,etax,
                    etam1,etam2,etam3,etam4,etam5,delta1,
                    tseq,B,seed){

  NIE_ss=matrix(data=0,nrow=5,ncol=length(tseq))
  S1_list=matrix(data=0,nrow=5,ncol=length(tseq)*32)
  S0_list=matrix(data=0,nrow=5,ncol=length(tseq)*32)
  
  for (ss in 1:5){
    
    
    ###############Generate S0/S1 for preparing the NIE_sub###############
    ###S1/S0 for Z of nth event is T1
    S_1<-NIE_prep(X,X2,n,T0,T1,sigmav,tau,
                  lambda.m1,lambda.m2,lambda.m3,lambda.m4,lambda.m5,
                  t.jump.m1,t.jump.m2,t.jump.m3,t.jump.m4,t.jump.m5,
                  lambda.y,t.jump.y,
                  betaz1,betaz2,betaz3,betaz4,betaz5,
                  betax1,betax2,betax3,betax4,betax5,
                  delta_X1,delta_X2,delta_X3,delta_X4,delta_X5,
                  etaz,etax,
                  etam1,etam2,etam3,etam4,etam5,delta1,
                  ss,T1,tseq,B,seed)[[2]]
    ###S1/S0 for Z of nth event is T0
    S_0<-NIE_prep(X,X2,n,T0,T1,sigmav,tau,
                  lambda.m1,lambda.m2,lambda.m3,lambda.m4,lambda.m5,
                  t.jump.m1,t.jump.m2,t.jump.m3,t.jump.m4,t.jump.m5,
                  lambda.y,t.jump.y,
                  betaz1,betaz2,betaz3,betaz4,betaz5,
                  betax1,betax2,betax3,betax4,betax5,
                  delta_X1,delta_X2,delta_X3,delta_X4,delta_X5,
                  etaz,etax,
                  etam1,etam2,etam3,etam4,etam5,delta1,
                  ss,T0,tseq,B,seed)[[2]]
    if (ss==1){
      S000000_n=S_0[1:length(tseq)]
      S011111_n=S_1[(length(tseq)*16+1):(length(tseq)*17)]
      S100000_n=S_0[(length(tseq)*16+1):(length(tseq)*17)]
      S111111_n=S_1[(length(tseq)*31+1):(length(tseq)*32)]
    }
    #Survival trt is T1
    NIE10_1=S_0[(length(tseq)*16+1):(length(tseq)*17)]- S_1[(length(tseq)*16+1):(length(tseq)*17)]
    NIE10_2=S_0[(length(tseq)*17+1):(length(tseq)*18)]- S_1[(length(tseq)*17+1):(length(tseq)*18)]
    NIE10_3=S_0[(length(tseq)*18+1):(length(tseq)*19)]- S_1[(length(tseq)*18+1):(length(tseq)*19)]
    NIE10_4=S_0[(length(tseq)*19+1):(length(tseq)*20)]- S_1[(length(tseq)*19+1):(length(tseq)*20)]
    NIE10_5=S_0[(length(tseq)*20+1):(length(tseq)*21)]- S_1[(length(tseq)*20+1):(length(tseq)*21)]
    NIE10_6=S_0[(length(tseq)*21+1):(length(tseq)*22)]- S_1[(length(tseq)*21+1):(length(tseq)*22)]
    NIE10_7=S_0[(length(tseq)*22+1):(length(tseq)*23)]- S_1[(length(tseq)*22+1):(length(tseq)*23)]
    NIE10_8=S_0[(length(tseq)*23+1):(length(tseq)*24)]- S_1[(length(tseq)*23+1):(length(tseq)*24)]
    NIE10_9=S_0[(length(tseq)*24+1):(length(tseq)*25)]- S_1[(length(tseq)*24+1):(length(tseq)*25)]
    NIE10_10=S_0[(length(tseq)*25+1):(length(tseq)*26)]- S_1[(length(tseq)*25+1):(length(tseq)*26)]
    NIE10_11=S_0[(length(tseq)*26+1):(length(tseq)*27)]- S_1[(length(tseq)*26+1):(length(tseq)*27)]
    NIE10_12=S_0[(length(tseq)*27+1):(length(tseq)*28)]- S_1[(length(tseq)*27+1):(length(tseq)*28)]
    NIE10_13=S_0[(length(tseq)*28+1):(length(tseq)*29)]- S_1[(length(tseq)*28+1):(length(tseq)*29)]
    NIE10_14=S_0[(length(tseq)*29+1):(length(tseq)*30)]- S_1[(length(tseq)*29+1):(length(tseq)*30)]
    NIE10_15=S_0[(length(tseq)*30+1):(length(tseq)*31)]- S_1[(length(tseq)*30+1):(length(tseq)*31)]
    NIE10_16=S_0[(length(tseq)*31+1):(length(tseq)*32)]- S_1[(length(tseq)*31+1):(length(tseq)*32)]
    
    NIE10=(NIE10_1+NIE10_2+NIE10_3+NIE10_4+NIE10_5+NIE10_6+NIE10_7+NIE10_8+NIE10_9
              +NIE10_10+NIE10_11+NIE10_12+NIE10_13+NIE10_14+NIE10_15+NIE10_16)/16
    #Survival trt is T0
    NIE01_1=S_1[(1:length(tseq))]- S_0[(1:length(tseq))]
    NIE01_2=S_1[(length(tseq)+1):(length(tseq)*2)]- S_0[(length(tseq)+1):(length(tseq)*2)]
    NIE01_3=S_1[(length(tseq)*2+1):(length(tseq)*3)]- S_0[(length(tseq)*2+1):(length(tseq)*3)]
    NIE01_4=S_1[(length(tseq)*3+1):(length(tseq)*4)]- S_0[(length(tseq)*3+1):(length(tseq)*4)]
    NIE01_5=S_1[(length(tseq)*4+1):(length(tseq)*5)]- S_0[(length(tseq)*4+1):(length(tseq)*5)]
    NIE01_6=S_1[(length(tseq)*5+1):(length(tseq)*6)]- S_0[(length(tseq)*5+1):(length(tseq)*6)]
    NIE01_7=S_1[(length(tseq)*6+1):(length(tseq)*7)]- S_0[(length(tseq)*6+1):(length(tseq)*7)]
    NIE01_8=S_1[(length(tseq)*7+1):(length(tseq)*8)]- S_0[(length(tseq)*7+1):(length(tseq)*8)]
    NIE01_9=S_1[(length(tseq)*8+1):(length(tseq)*9)]- S_0[(length(tseq)*8+1):(length(tseq)*9)]
    NIE01_10=S_1[(length(tseq)*9+1):(length(tseq)*10)]- S_0[(length(tseq)*9+1):(length(tseq)*10)]
    NIE01_11=S_1[(length(tseq)*10+1):(length(tseq)*11)]- S_0[(length(tseq)*10+1):(length(tseq)*11)]
    NIE01_12=S_1[(length(tseq)*11+1):(length(tseq)*12)]- S_0[(length(tseq)*11+1):(length(tseq)*12)]
    NIE01_13=S_1[(length(tseq)*12+1):(length(tseq)*13)]- S_0[(length(tseq)*12+1):(length(tseq)*13)]
    NIE01_14=S_1[(length(tseq)*13+1):(length(tseq)*14)]- S_0[(length(tseq)*13+1):(length(tseq)*14)]
    NIE01_15=S_1[(length(tseq)*14+1):(length(tseq)*15)]- S_0[(length(tseq)*14+1):(length(tseq)*15)]
    NIE01_16=S_1[(length(tseq)*15+1):(length(tseq)*16)]- S_0[(length(tseq)*15+1):(length(tseq)*16)]
    
    NIE01=(NIE01_1+NIE01_2+NIE01_3+NIE01_4+NIE01_5+NIE01_6+NIE01_7+NIE01_8+NIE01_9
              +NIE01_10+NIE01_11+NIE01_12+NIE01_13+NIE01_14+NIE01_15+NIE01_16)/16
    
    NIE_ss[ss,]=(NIE01-NIE10)/2
    
    S1_list[ss,]<-S_1
    S0_list[ss,]<-S_0
  }
  
NDE01=S100000_n-S000000_n
NDE10=S011111_n-S111111_n
NIE01=S011111_n-S000000_n
NIE10=S100000_n-S111111_n

NDE=(NDE01-NDE10)/2
NIE=(NIE01-NIE10)/2
TE=NDE+NIE

NIE1<-NIE_ss[1,]
NIE2<-NIE_ss[2,]
NIE3<-NIE_ss[3,]
NIE4<-NIE_ss[4,]
NIE5<-NIE_ss[5,]

effect<-as.data.frame(list(NDE=NDE,NIE=NIE,TE=TE,NIE1=NIE1,NIE2=NIE2,NIE3=NIE3,
                           NIE4=NIE4,NIE5=NIE5))
return(list(effect=effect,S1=S1_list, S0=S0_list))

}

