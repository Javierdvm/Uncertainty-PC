#Limpiar
rm(list=ls())
#Paquetes
library(MASS)
library(matrixStats)
library(Matrix)
#Parámetros iniciales
R=1
B=500
t=100
N=100
r=1
phi1<-0.5
ro<-0
Tau<-1
significance=0.05
z<-qnorm(1-significance/2,0,1)
Fi<-matrix(0,r,t)
PC1ESTINV<-matrix(0,t,r)
PC1EST<-matrix(0,t,r)
F<-matrix(0,r,t)
FA<-matrix(0,r,t)
Dif1<-matrix(0,1,r)
Dif2<-matrix(0,1,r)
PC1BAI<-matrix(0,t,r)
Load<-matrix(0,N,r)
E<-matrix(0,N,t)
A<-matrix(0,N,t)
LLL_ALL<-matrix(0,r,r)
cover=matrix(0,t,1)
varPC1<-matrix(0,t,1)
varPC2<-matrix(0,t,1)
varPC3<-matrix(0,t,1)
CIL<-matrix(0,t,1)
CIU<-matrix(0,t,1)
var_asin<-matrix(0,t,r)
desv_asin<-matrix(0,t,r)
Error_asin<-matrix(0,t,r)
MError_asin<-matrix(0,R,r)
Ft_boot<-matrix(0,t,B)
tPC1BAI1_boot<-matrix(0,t,B)
MError_boot<-matrix(0,r,R)
MSE<-matrix(0,t,r)
RMSE<-matrix(0,R,r)
tPC1boot<-matrix(0,t,B)
RMSEB<-matrix(0,R,r)
Ft_b<-matrix(0,t,r)
ar<-matrix(0,B,r)
ar_boot<-matrix(0,B,r)

COVERAGES<-rep(0,t)
P<-matrix(runif(r*N,0,1),nrow=N,ncol=r)
PI<-ro*diag(N)

for (jj in 1:R) {
  print(jj)
  #Crear Factor
  Fi[,1]<-matrix(rnorm(r,0,(1/1-phi1^2)),r,1)
  E[,1]<-matrix(rnorm(N,0,(1/1-ro^2)),N,1)
  eta<-matrix(rnorm(r,0,1),r,1)
  for(ii in 2:t) {
    eta<-matrix(rnorm(r,0,1),r,1)
    Fi[,ii]<-matrix(phi1%*%Fi[,ii-1]+sqrt((1-phi1^2))*eta,r,1)
    #A[,ii]<-matrix(rnorm(N,0,Tau),nrow=N,ncol=1)
    A[,ii]<-matrix(mvrnorm(1,rep(0,N),Tau*diag(N)),nrow=N,ncol=1)
    E[,ii]<-PI%*%E[,ii-1]+(diag(N)-PI^2)^(1/2)%*%A[,ii]
  }
  F<-matrix(Fi,nrow=r,ncol=t)
  #Crear Variables e idiosincrático
  y<-matrix(P%*%F+E,nrow=N,ncol=t)
  
  #BAI
  Y<-t(y)
  RR<-Y%*%t(Y)
  eR<-eigen(RR)
  values<-eR$values[c(1:r)]
  vectors<-matrix(eR$vectors[,c(1:r)],t,r)
  PC1EST<-sqrt(t)*vectors
  #Efecto espejo
  PC1ESTINV<-(-1)*PC1EST
  Dif1<-colSums(abs(t(F)-PC1EST))
  Dif2<-colSums(abs(t(F)-PC1ESTINV))
  if (Dif1<Dif2) {
    PC1R=PC1EST
  } else {
    PC1R=PC1ESTINV
  }
  tPC1<-ts(PC1R)
  Load<-matrix(t(Y)%*%PC1R/t,N,r)
  Resid<-y-Load%*%t(PC1R)
  
  #OLS
  fit<-ar.ols(ts(tPC1),aic=FALSE,order.max=1,demean=FALSE,intercept=FALSE)
  resid<-scale(fit$resid,center=TRUE,scale=FALSE)
  resid[is.na(resid)]<-mean(resid[-1])

  #Bootstrap
  for (bb in 1:B) {
    resid_boot<-sample(resid,size = t,replace = TRUE,prob = NULL)
    Ft_boot[1,]=tPC1[1]
    Ft_b[1]=tPC1[1]
    for (rr in 2:t) {
     Ft_b[rr]<-fit$ar*Ft_b[rr-1]+resid_boot[rr]
    }
    fit2<-ar.ols(ts(Ft_b),aic=FALSE,order.max=1,demean=FALSE,intercept=FALSE)
    ar[bb]<-fit2$ar
  }
  ar_boot<-sample(ar,size = B,replace = TRUE,prob = NULL)
  for (bb in 1:B) {
    resid_boot<-sample(resid,size = t,replace = TRUE,prob = NULL)
    Ft_boot[1,]=tPC1[1]
    for (tt in 2:t) {
      Ft_boot[tt,bb]<-ar_boot[bb]*tPC1[tt-1]+resid_boot[tt]
    }
    tPC1BAI1_boot[,bb]<-ts(Ft_boot[,bb])
  
    plot(ts(Ft_boot[,bb]))
    lines(tPC1,col="red")
    lines(ts(t(F)),col="blue")
    
  }
  #Confidence Intervals
  #Varianza Bootstrap
  var_boot<-rowVars(tPC1BAI1_boot)
  desv_boot<-sqrt(var_boot)
  Error_boot<-desv_boot
  for (qq in 1:t-1) {
    CIL[qq]<-quantile(tPC1BAI1_boot[qq+1,],0.05)    
    CIU[qq]<-quantile(tPC1BAI1_boot[qq+1,],0.95)    
  }
  #CIL<-tPC1-z*Error_boot
  #CIU<-tPC1+z*Error_boot
  
  #RMSE
  MSE<-((t(F)-PC1R)^2)
  RMSE[jj]=sqrt(colMeans(MSE))
  var_boot<-rowVars(tPC1BAI1_boot)
  RMSEB[jj]<-sqrt(mean(var_boot))
  
  
  #Coverages
  cover[1] <- if (CIL[1]<F[1]&CIU[1]>F[1]) cover[1]=0 else cover[1]=1
  for(zz in 2:t-1) {
    cover[zz] <- if (CIL[zz]<F[zz]&CIU[zz]>F[zz]) cover[zz]=0 else cover[zz]=1
  }
  COVER<-matrix(cover,nrow=t,ncol=r)
  COVERAGES=COVERAGES+COVER
  
  #Gráficos Factor
  plot(ts(t(F)))
  lines(tPC1,col="blue")
  lines(ts(CIL),col="red")
  lines(ts(CIU),col="red")
  
}
print(t)
print(N)
print(phi1)
print(Tau)
print(ro)
MRMSE<-colMeans(RMSE)
MRMSEB<-colMeans(RMSEB)
print(MRMSE)
print(MRMSEB)
coverage_ratio=1-colSums(COVERAGES)/(t*R)
print(coverage_ratio)
print(fit$ar)
