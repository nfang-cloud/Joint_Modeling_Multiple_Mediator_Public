
####data generating with simulated X,Z


datagen_X<-function(n,vi,tau,betaz1,betaz2,betax1,betax2,lambda.m1,lambda.m2,t.jump.m1,t.jump.m2,
                    etaz,etax,etam1,etam2,delta1,deltax1,deltax2,delta2_1,delta2_2,
                    lambda.y,t.jump.y,cen=TRUE){
  
  #Generate X and Z
  #Generate the corrlation matrix
  k=NROW(betax1)
  mu<- rep(0,k)
  corMat <- rcorrmatrix(k)
  
  # Generate data from normal distribution
  X <- mvrnorm(n, mu=mu, Sigma=corMat, empirical = TRUE)
  
  # Generate binary trt
  Z <- rbinom(n,1,0.5)
  p=ncol(X)
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all1=c(betaz1,betax1,deltax1)
  risk1_0=c(exp(cbind(0,X,vi)%*%beta.all1))
  risk1_1=c(exp(cbind(1,X,vi)%*%beta.all1))
  
  beta.all2=c(betaz2,betax2,deltax2)
  risk2_0=c(exp(cbind(0,X,vi)%*%beta.all2))
  risk2_1=c(exp(cbind(1,X,vi)%*%beta.all2))
  ####generate survival time until all individual has follow-up at lest time tau
  M1_0=mysimRec(t.jump.m1,lambda.m1,tau,risk1_0,recurrent=TRUE)
  M1_1=mysimRec(t.jump.m1,lambda.m1,tau,risk1_1,recurrent=TRUE)
  
  M2_0=mysimRec(t.jump.m2,lambda.m2,tau,risk2_0,recurrent=TRUE)
  M2_1=mysimRec(t.jump.m2,lambda.m2,tau,risk2_1,recurrent=TRUE)
  
  M1_0$event_Tp=1
  M1_1$event_Tp=1
  M2_0$event_Tp=2
  M2_1$event_Tp=2
  M1_0[which(M1_0$event==0),]$event_Tp=0
  M1_1[which(M1_1$event==0),]$event_Tp=0
  M2_0[which(M2_0$event==0),]$event_Tp=0
  M2_1[which(M2_1$event==0),]$event_Tp=0
  
  M_00=rbind(M1_0,M2_0)
  M_01=rbind(M1_0,M2_1)
  M_10=rbind(M1_1,M2_0)
  M_11=rbind(M1_1,M2_1)
  
  M_00=M_00[order(M_00$ID,M_00$stoptime),]
  M_01=M_01[order(M_01$ID,M_01$stoptime),]
  M_10=M_10[order(M_10$ID,M_10$stoptime),]
  M_11=M_11[order(M_11$ID,M_11$stoptime),]
  
  M_00=M_00[!duplicated( M_00[ , c("ID","stoptime","event","event_Tp")]),]
  M_01=M_01[!duplicated( M_01[ , c("ID","stoptime","event","event_Tp")]),]
  M_10=M_10[!duplicated( M_10[ , c("ID","stoptime","event","event_Tp")]),]
  M_11=M_11[!duplicated( M_11[ , c("ID","stoptime","event","event_Tp")]),]
  
  
  ####Simulate T^{zM^z'} one at a time
  etam1=c(etam1,delta2_1)
  etam2=c(etam2,delta2_2)
  T000=T001=T010=T011=T100=T101=T110=T111=rep(0,n)
  risk0.y=c(exp(cbind(0,X,vi)%*%c(etaz,etax,delta1)))
  risk1.y=c(exp(cbind(1,X,vi)%*%c(etaz,etax,delta1)))
  
  for (i in 1:n){
    tmpM_00_1=M_00[which(M_00$ID==i & M_00$event_Tp==1),]
    tmpM_00_2=M_00[which(M_00$ID==i & M_00$event_Tp==2),]
    tmpM_01_1=M_01[which(M_01$ID==i & M_01$event_Tp==1),]
    tmpM_01_2=M_01[which(M_01$ID==i & M_01$event_Tp==2),]
    tmpM_10_1=M_10[which(M_10$ID==i & M_10$event_Tp==1),]
    tmpM_10_2=M_10[which(M_10$ID==i & M_10$event_Tp==2),]
    tmpM_11_1=M_11[which(M_11$ID==i & M_11$event_Tp==1),]
    tmpM_11_2=M_11[which(M_11$ID==i & M_11$event_Tp==2),]
    
    
    lambda.risk00_1=exp((0:length(tmpM_00_1$stoptime))*(etam1[1]+etam1[2]*vi[i]))
    t.jump.risk00_1=tmpM_00_1$stoptime
    
    lambda.risk00_2=exp((0:length(tmpM_00_2$stoptime))*(etam2[1]+etam2[2]*vi[i]))
    t.jump.risk00_2=tmpM_00_2$stoptime
    
    lambda.risk01_1=exp((0:length(tmpM_01_1$stoptime))*(etam1[1]+etam1[2]*vi[i]))
    t.jump.risk01_1=tmpM_00_1$stoptime
    
    lambda.risk01_2=exp((0:length(tmpM_01_2$stoptime))*(etam2[1]+etam2[2]*vi[i]))
    t.jump.risk01_2=tmpM_01_2$stoptime
    
    lambda.risk10_1=exp((0:length(tmpM_10_1$stoptime))*(etam1[1]+etam1[2]*vi[i]))
    t.jump.risk10_1=tmpM_10_1$stoptime
    
    lambda.risk10_2=exp((0:length(tmpM_10_2$stoptime))*(etam2[1]+etam2[2]*vi[i]))
    t.jump.risk10_2=tmpM_10_2$stoptime
    
    lambda.risk11_1=exp((0:length(tmpM_11_1$stoptime))*(etam1[1]+etam1[2]*vi[i]))
    t.jump.risk11_1=tmpM_11_1$stoptime
    
    lambda.risk11_2=exp((0:length(tmpM_11_2$stoptime))*(etam2[1]+etam2[2]*vi[i]))
    t.jump.risk11_2=tmpM_11_2$stoptime
    
    lambda.risk00<-myprod(t.jump.risk00_1,lambda.risk00_1,t.jump.risk00_2,lambda.risk00_2)[[2]]
    t.jump.risk00<-myprod(t.jump.risk00_1,lambda.risk00_1,t.jump.risk00_2,lambda.risk00_2)[[1]]
    
    lambda.risk01<-myprod(t.jump.risk01_1,lambda.risk01_1,t.jump.risk01_2,lambda.risk01_2)[[2]]
    t.jump.risk01<-myprod(t.jump.risk01_1,lambda.risk01_1,t.jump.risk01_2,lambda.risk01_2)[[1]]
    
    lambda.risk10<-myprod(t.jump.risk10_1,lambda.risk10_1,t.jump.risk10_2,lambda.risk10_2)[[2]]
    t.jump.risk10<-myprod(t.jump.risk10_1,lambda.risk10_1,t.jump.risk10_2,lambda.risk10_2)[[1]]
    
    lambda.risk11<-myprod(t.jump.risk11_1,lambda.risk11_1,t.jump.risk11_2,lambda.risk11_2)[[2]]
    t.jump.risk11<-myprod(t.jump.risk11_1,lambda.risk11_1,t.jump.risk11_2,lambda.risk11_2)[[1]]
    
    my00=myprod(t.jump.y,lambda.y,t.jump.risk00,lambda.risk00)
    my01=myprod(t.jump.y,lambda.y,t.jump.risk01,lambda.risk01)
    my10=myprod(t.jump.y,lambda.y,t.jump.risk10,lambda.risk10)
    my11=myprod(t.jump.y,lambda.y,t.jump.risk11,lambda.risk11)
    
    ###bring the notation 0/1 to front 
    T000[i]=mysimRec(my00[[1]],my00[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T100[i]=mysimRec(my00[[1]],my00[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
    
    T010[i]=mysimRec(my10[[1]],my10[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T110[i]=mysimRec(my10[[1]],my10[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
    
    T001[i]=mysimRec(my01[[1]],my01[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T101[i]=mysimRec(my01[[1]],my01[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
    
    T011[i]=mysimRec(my11[[1]],my11[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T111[i]=mysimRec(my11[[1]],my11[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
  }
  ####Simulate Z
  ID1=ID[which(Z==1)]
  ID0=ID[which(Z==0)]
  Mdat=rbind(M_00[which(M_00$ID%in%ID0),],M_11[which(M_11$ID%in%ID1),])
  Mdat=Mdat[which(Mdat$stoptime<=tau),]
  Mdat=Mdat[order(Mdat$ID,Mdat$stoptime),]
  T=T000*(1-Z)+T111*Z
  C=tau-0.01
  if (cen){
    C=runif(n,tau/2,tau)
  }
  Y=pmin(C,T)
  D=as.numeric(T<=C)
  Tp=rep(0,n)
  ydata=data.frame(ID,Y,D,Tp,Z,X)
  names(ydata)=c("ID","Y","D","Tp",paste("X",as.character(1:(p+1)),sep=""))
  mdata=Mdat
  adata=merge(ydata,mdata,by="ID")
  adata=adata[which(adata$stoptime<adata$Y),]
  adata$D=adata$event
  adata$Y=adata$stoptime
  adata$Tp=adata$event_Tp
  adata=adata[,c("ID","Y","D","Tp",paste("X",as.character(1:(p+1)),sep=""))]
  ydata$D=ydata$D*2
  alldata=rbind(adata,ydata)
  alldata=alldata[order(alldata$ID,alldata$Y),]
  ###first is Z, all others are X
  names(alldata)=c("ID","stoptime","event","event_Tp",paste("X",as.character(1:(p+1)),sep=""))
  ydata$vi=vi
  data=list(alldata=alldata,edat=data.frame(ID,T000,T100,T001,T101,T010,T110,T011,T111),mdat=merge(data.frame(Mdat),ydata[,-c(2,3)],by="ID"))

  }
