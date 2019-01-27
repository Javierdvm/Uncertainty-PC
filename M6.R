#Limpiar
rm(list=ls())
#Paquetes
library(MASS)
library(matrixStats)
library(Matrix)
#Parámetros iniciales
t=150
N=200
ro<-0

B=1000
R=100
r=1

significance1=0.2
significance2=0.05
significance3=0.01
significance4=0.2
significance5=0.05
significance6=0.01
z1<-qnorm(1-significance1/2,0,1)
z2<-qnorm(1-significance2/2,0,1)
z3<-qnorm(1-significance3/2,0,1)
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
LLL_ALL<-matrix(0,r,r)
cover1=matrix(0,t,1)
cover2=matrix(0,t,1)
cover3=matrix(0,t,1)
cover4=matrix(0,t,1)
cover5=matrix(0,t,1)
cover6=matrix(0,t,1)
varPC1<-matrix(0,t,1)
varPC2<-matrix(0,t,1)
varPC3<-matrix(0,t,1)
CIL1<-matrix(0,t,1)
CIU1<-matrix(0,t,1)
CIL2<-matrix(0,t,1)
CIU2<-matrix(0,t,1)
CIL3<-matrix(0,t,1)
CIU3<-matrix(0,t,1)
CIL4<-matrix(0,t,1)
CIU4<-matrix(0,t,1)
CIL5<-matrix(0,t,1)
CIU5<-matrix(0,t,1)
CIL6<-matrix(0,t,1)
CIU6<-matrix(0,t,1)
differences<-matrix(0,t,r)
Error_asin<-matrix(0,t,r)
MSE<-matrix(0,t,r)
RMSE<-matrix(0,R,r)
RMSEB<-matrix(0,R,r)
difRMS<-matrix(0,R,r)
difff<-matrix(0,12,1)
MRMSE<-matrix(0,12,1)
MRMSEB<-matrix(0,12,1)
coverage_ratio1<-matrix(0,3,1)
coverage_ratio2<-matrix(0,3,1)
coverage_ratio3<-matrix(0,3,1)
coverage_ratio4<-matrix(0,3,1)
coverage_ratio5<-matrix(0,3,1)
coverage_ratio6<-matrix(0,3,1)
q<-matrix(0,12,1)
phi<-matrix(0,12,1)
A<-matrix(0,N,t)
var_asin<-matrix(0,t,r)
desv_asin<-matrix(0,t,r)
Error_asin<-matrix(0,t,r)
MError_asin<-matrix(0,R,r)
Ft_boot<-matrix(0,t,B)
tPC1BAI1_boot<-matrix(0,t,B)
MError_boot<-matrix(0,r,R)
tPC1boot<-matrix(0,t,B)
RMSEB<-matrix(0,R,r)
Ft_b<-matrix(0,t,r)
ar<-matrix(0,B,r)
ar_boot<-matrix(0,B,r)
arr<-matrix(0,1,R)


rrr=0
for (SN in c(10,1,0.5,0.1)){
  Tau<-1/SN
  for (phi1 in c(0.2,0.5,0.9)) {
    rrr=rrr+1
    COVERAGES1<-rep(0,t)
    COVERAGES2<-rep(0,t)
    COVERAGES3<-rep(0,t)
    COVERAGES4<-rep(0,t)
    COVERAGES5<-rep(0,t)
    COVERAGES6<-rep(0,t)
    P<-matrix(runif(r*N,0,1),nrow=N,ncol=r)
    PI<-ro*diag(N)
    U<-matrix(diag(runif(N,0.8,1.2)),nrow=N, ncol=N)
    for (jj in 1:R) {
      print(SN)
      print(phi1)
      print(jj)

#Crear Factor
Fi[,1]<-matrix(rnorm(r,0,(1/1-phi1^2)),r,1)
E[,1]<-matrix(rnorm(N,0,(1/1-ro^2)),N,1)
eta<-matrix(rnorm(r,0,1),r,1)
for(ii in 2:t) {
  eta<-matrix(rnorm(r,0,1),r,1)
  Fi[,ii]<-matrix(phi1%*%Fi[,ii-1]+sqrt((1-phi1^2))*eta,r,1)
  A[,ii]<-matrix(mvrnorm(1,rep(0,N),Tau*diag(N)),nrow=N,ncol=1)
 E[,ii]<-PI%*%E[,ii-1]+(diag(N)-PI^2)^(1/2)%*%A[,ii]
}
F<-matrix(Fi,nrow=r,ncol=t)
#Crear Variables e idiosincrático
expo<-1:N
Parametro<-0.5*(Tau^2)
Psi<-head((0.5/expo)*(Tau^2),N)
Psi<-replace(Psi, Psi==0.5*(Tau^2), Tau)
covar<-toeplitz(Psi)
#E<-matrix(mvrnorm(t,rep(0,N),covar),nrow=N,ncol=t)
#E<-matrix(mvrnorm(t,rep(0,N),Tau*U),nrow=N,ncol=t)
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
  PC1=PC1EST
} else {
  PC1=PC1ESTINV
}
tPC1BAI1<-ts(PC1)
Load<-matrix(t(Y)%*%PC1/t,N,r)
Resid<-y-Load%*%t(PC1)
tPC1<-ts(PC1)

  #OLS
  fit<-ar.ols(ts(PC1),aic=FALSE,order.max=1,demean=FALSE,intercept=FALSE)
  resid<-scale(fit$resid,center=FALSE,scale=FALSE)
  resid[is.na(resid)]<-mean(resid[-1])
  arr[jj]<-fit$ar
  #Bootstrap AR
  for (bb in 1:B) {
    resid_boot<-sample(resid,size = t,replace = TRUE,prob = NULL)
    Ft_boot[1,]=tPC1[1]
    Ft_b[1]=ts(PC1)[1]
    for (rr in 2:t) {
     Ft_b[rr]<-fit$ar*Ft_b[rr-1]+resid_boot[rr]
    }
    fit2<-ar.ols(ts(Ft_b),aic=FALSE,order.max=1,demean=FALSE,intercept=FALSE)
    ar[bb]<-fit2$ar
  }
 Resid<-scale(Resid,center=FALSE,scale=FALSE)
#Bootstrap Idiosincrático
  for (bb in 1:B) {
    resamplet <- sample(nrow(Resid),size = N,replace = TRUE,prob = NULL)
    Resid_boot<-Resid[resamplet,]
    for (tt in 2:t) {
      Ft_boot[tt,bb]<-ar[bb]*tPC1[tt-1]+resid[tt]
    }
    y_boot<-matrix(Load%*%t(Ft_boot[,bb])+Resid_boot,nrow=N,ncol=t)
    
#Sacando PC del bootstrap
    #BAI
    Y<-t(y_boot)
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
      PC1B=PC1EST
    } else {
      PC1B=PC1ESTINV
    }
    tPC1BAI1_boot[,bb]<-ts(PC1B)
    Load_boot<-matrix(t(Y)%*%PC1B/t,N,r)
  
  }
  
  #Varianza Bootstrap, mse
  var_boot<-rowVars(tPC1BAI1_boot)
  desv_boot<-sqrt(var_boot)
  Error_boot<-desv_boot
  RMSEB[jj]<-sqrt(mean(var_boot))
for (qq in 1:t) {
  differences[qq]<-(var_boot[qq])-((F[qq]-PC1[qq])^2)
}
MSE<-((t(F)-PC1)^2)
RMSE[jj]=sqrt(colMeans(MSE))
difRMS[jj]<-colMeans(differences)

#Confidence Intervals
for (qq in 1:t) {
#IC1
CIL1[qq]<-PC1[qq]-z1*Error_boot[qq]
CIU1[qq]<-PC1[qq]+z1*Error_boot[qq]
#IC2
CIL2[qq]<-PC1[qq]-z2*Error_boot[qq]
CIU2[qq]<-PC1[qq]+z2*Error_boot[qq]
#IC3
CIL3[qq]<-PC1[qq]-z3*Error_boot[qq]
CIU3[qq]<-PC1[qq]+z3*Error_boot[qq]
} 
#CI4
for (qq in 1:t) {
  CIL4[qq]<-quantile(tPC1BAI1_boot[qq,],significance4/2)    
  CIU4[qq]<-quantile(tPC1BAI1_boot[qq,],1-significance4/2)    
}
#CI5
for (qq in 1:t) {
  CIL5[qq]<-quantile(tPC1BAI1_boot[qq,],significance5/2)    
  CIU5[qq]<-quantile(tPC1BAI1_boot[qq,],1-significance5/2)    
}
#CI6
for (qq in 1:t) {
  CIL6[qq]<-quantile(tPC1BAI1_boot[qq,],significance6/2)    
  CIU6[qq]<-quantile(tPC1BAI1_boot[qq,],1-significance6/2)    
}


#Coverages1
cover1[1] <- if (CIL1[1]<F[1]&CIU1[1]>F[1]) cover1[1]=0 else cover1[1]=1
for(zz in 2:t) {
  cover1[zz] <- if (CIL1[zz]<F[zz]&CIU1[zz]>F[zz]) cover1[zz]=0 else cover1[zz]=1
}
COVER1<-matrix(cover1,nrow=t,ncol=r)
COVERAGES1=COVERAGES1+COVER1

#coverages2
cover2[1] <- if (CIL2[1]<F[1]&CIU2[1]>F[1]) cover2[1]=0 else cover2[1]=1
for(zz in 2:t) {
  cover2[zz] <- if (CIL2[zz]<F[zz]&CIU2[zz]>F[zz]) cover2[zz]=0 else cover2[zz]=1
}
COVER2<-matrix(cover2,nrow=t,ncol=r)
COVERAGES2=COVERAGES2+COVER2
#coverages3
cover3[1] <- if (CIL3[1]<F[1]&CIU3[1]>F[1]) cover3[1]=0 else cover3[1]=1
for(zz in 2:t) {
  cover3[zz] <- if (CIL3[zz]<F[zz]&CIU3[zz]>F[zz]) cover3[zz]=0 else cover3[zz]=1
}
COVER3<-matrix(cover3,nrow=t,ncol=r)
COVERAGES3=COVERAGES3+COVER3
#coverages4
cover4[1] <- if (CIL4[1]<F[1]&CIU4[1]>F[1]) cover4[1]=0 else cover4[1]=1
for(zz in 2:t) {
  cover4[zz] <- if (CIL4[zz]<F[zz]&CIU4[zz]>F[zz]) cover4[zz]=0 else cover4[zz]=1
}
COVER4<-matrix(cover4,nrow=t,ncol=r)
COVERAGES4=COVERAGES4+COVER4
#coverages5
cover5[1] <- if (CIL5[1]<F[1]&CIU5[1]>F[1]) cover5[1]=0 else cover5[1]=1
for(zz in 2:t) {
  cover5[zz] <- if (CIL5[zz]<F[zz]&CIU5[zz]>F[zz]) cover5[zz]=0 else cover5[zz]=1
}
COVER5<-matrix(cover5,nrow=t,ncol=r)
COVERAGES5=COVERAGES5+COVER5
#coverages6
cover6[1] <- if (CIL6[1]<F[1]&CIU6[1]>F[1]) cover6[1]=0 else cover6[1]=1
for(zz in 2:t) {
  cover6[zz] <- if (CIL6[zz]<F[zz]&CIU6[zz]>F[zz]) cover6[zz]=0 else cover6[zz]=1
}
COVER6<-matrix(cover6,nrow=t,ncol=r)
COVERAGES6=COVERAGES6+COVER6
#Gráficos Factor
plot(ts(t(F)))
lines(tPC1BAI1,col="blue")
lines(ts(CIL1),col="red")
lines(ts(CIU1),col="red")


}
MRMSE[rrr]<-colMeans(RMSE)
MRMSEB[rrr]<-colMeans(RMSEB)
difff[rrr]<-colMeans(difRMS)
coverage_ratio1[rrr]=1-colSums(COVERAGES1)/(t*R)
coverage_ratio2[rrr]=1-colSums(COVERAGES2)/(t*R)
coverage_ratio3[rrr]=1-colSums(COVERAGES3)/(t*R)
coverage_ratio4[rrr]=1-colSums(COVERAGES4)/(t*R)
coverage_ratio5[rrr]=1-colSums(COVERAGES5)/(t*R)
coverage_ratio6[rrr]=1-colSums(COVERAGES6)/(t*R)
q[rrr]<-SN
phi[rrr]<-phi1
}
}
#hist(ar)
#print(fit$ar)
#print(mean(ar))
salida <- data.frame(q,phi,MRMSE,MRMSEB,difff,coverage_ratio1, coverage_ratio2,coverage_ratio3,coverage_ratio4,coverage_ratio5,coverage_ratio6)
write.csv2(salida,"M6 b=1000 T=150 N=200.csv")
