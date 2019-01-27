#Limpiar
rm(list=ls())
#Paquetes
library(MASS)
library(matrixStats)
library(Matrix)
#Parámetros iniciales
t=100
N=200
ro<-0

B=75
R=20
r=2

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
RMSE1<-matrix(0,R,1)
RMSEB1<-matrix(0,R,1)
RMSE2<-matrix(0,R,1)
RMSEB2<-matrix(0,R,1)
difRMS<-matrix(0,R,r)
MRMSE2<-matrix(0,9,1)
MRMSE1<-matrix(0,9,1)
MRMSEB1<-matrix(0,9,1)
MRMSEB2<-matrix(0,9,1)
coverage_ratio1<-matrix(0,3,1)
coverage_ratio2<-matrix(0,3,1)
coverage_ratio3<-matrix(0,3,1)
coverage_ratio4<-matrix(0,3,1)
coverage_ratio5<-matrix(0,3,1)
coverage_ratio6<-matrix(0,3,1)
www<-matrix(0,9,1)
phi<-matrix(0,9,1)
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
tPC1b<-matrix(0,t,r)
ar<-matrix(0,B,r)
ar_boot<-matrix(0,B,r)
arr<-matrix(0,1,R)
Resid_boot2<-matrix(0,N,t)
Yt_b<-matrix(0,N,t)
Y_boot1<-matrix(0,t,N)
PC_boot1_inv<-matrix(0,t,r)
PC_boot1R<-matrix(0,t,r)
PC_boot2R<-matrix(0,t,r)
PC1_BOOT<-matrix(0,t,B)
PC2_BOOT<-matrix(0,t,B)
PC_pseudo2<-matrix(0,t,r)


rrr=0
for (SN in c(5,1,0.2)){ 
  Tau<-1/SN
  for (phi1 in c(0.2,0.5,0.8)) {
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
      for (q in 1:r) {
        F[q,1]<-rnorm(1,0,(1/1-phi1^2))
        for(ii in 2:t) {
          eta<-rnorm(1,0,1)
          F[q,ii]<-phi1%*%F[q,ii-1]+sqrt((1-phi1^2))*eta
        }
      }
      E[,1]<-matrix(rnorm(N,0,(1/1-ro^2)),N,1)
      for(ii in 2:t) {
        A[,ii]<-matrix(mvrnorm(1,rep(0,N),Tau*diag(N)),nrow=N,ncol=1)
        E[,ii]<-PI%*%E[,ii-1]+(diag(N)-PI^2)^(1/2)%*%A[,ii]
      }
      #Correlaciones y Heterocedasticidad
      expo<-1:N
      Parametro<-0.5*(Tau^2)
      Psi<-head((0.5/expo)*(Tau^2),N)
      Psi<-replace(Psi, Psi==0.5*(Tau^2), Tau)
      covar<-toeplitz(Psi)
      #E<-matrix(mvrnorm(t,rep(0,N),covar),nrow=N,ncol=t)
      #E<-matrix(mvrnorm(t,rep(0,N),Tau*U),nrow=N,ncol=t)
      y<-matrix(P%*%F+E,nrow=N,ncol=t)
      PCR<-matrix(0,t,r)
      #BAI
      Y<-t(y)
      RR<-Y%*%t(Y)
      eR<-eigen(RR)
      values<-eR$values[c(1:r)]
      sum(eR$values[c(1:r)])/sum(eR$values)
      vectors<-matrix(eR$vectors[,c(1:r)],t,r)
      PC<-sqrt(t)*vectors
      PC_INV1<-(-1)*PC
      Load<-matrix(t(Y)%*%PC/t,N,r)
      PCr<-PC%*%solve(solve(t(Load)%*%Load)%*%(t(Load)%*%P))
      Load<-matrix(t(Y)%*%PCr/t,N,r)
      Resid<-t(Y)-Load%*%t(PCr)
      PC_INV2<-(-1)*PCr
      
      if (r==1) {
        aa<-lm(F[1,] ~ 0+PCr[,1])
        bbb<-lm(F[1,] ~ 0+PC_INV2[,1])
        ee<-lm(F[1,] ~ 0+PC[,1])
        ff<-lm(F[1,] ~ 0+PC_INV1[,1])
    
      
        RS1<-matrix(c(summary(aa)$coefficients[1, 1],summary(bbb)$coefficients[1, 1],summary(ee)$coefficients[1, 1],summary(ff)$coefficients[1, 1]),4,1)
        
        AA = which.max(RS1)
        switch(AA,"1"={PCR[,1]=PCr[,1]},"2"={PCR[,1]=PC_INV2[,1]},"3"={PCR[,1]=PC[,1]},"4"={PCR[,1]=PC_INV1[,1]})
        
        
      } else {
      
      aa<-lm(F[1,] ~ 0+PCr[,1])
      bbb<-lm(F[1,] ~ 0+PC_INV2[,1])
      cc<-lm(F[1,] ~ 0+PCr[,2])
      dd<-lm(F[1,] ~ 0+PC_INV2[,2])
      ee<-lm(F[1,] ~ 0+PC[,1])
      ff<-lm(F[1,] ~ 0+PC_INV1[,1])
      gg<-lm(F[1,] ~ 0+PC[,2])
      hh<-lm(F[1,] ~ 0+PC_INV1[,2])
      
      ii<-lm(F[2,] ~ 0+PCr[,1])
      jjj<-lm(F[2,] ~ 0+PC_INV2[,1])
      kk<-lm(F[2,] ~ 0+PCr[,2])
      ll<-lm(F[2,] ~ 0+PC_INV2[,2])
      mm<-lm(F[2,] ~ 0+PC[,1])
      nn<-lm(F[2,] ~ 0+PC_INV1[,1])
      oo<-lm(F[2,] ~ 0+PC[,2])
      pp<-lm(F[2,] ~ 0+PC_INV1[,2])
      
      RS1<-matrix(c(summary(aa)$coefficients[1, 1],summary(bbb)$coefficients[1, 1],summary(cc)$coefficients[1, 1],summary(dd)$coefficients[1, 1],summary(ee)$coefficients[1, 1],summary(ff)$coefficients[1, 1],summary(gg)$coefficients[1, 1],summary(hh)$coefficients[1, 1]),8,1)
      RS2<-matrix(c(summary(ii)$coefficients[1, 1],summary(jjj)$coefficients[1, 1],summary(kk)$coefficients[1, 1],summary(ll)$coefficients[1, 1],summary(mm)$coefficients[1, 1],summary(nn)$coefficients[1, 1],summary(oo)$coefficients[1, 1],summary(pp)$coefficients[1, 1]),8,1)
      
      AA= which.max(RS1)
      BB= which.max(RS2)
      switch(AA,"1"={PCR[,1]=PCr[,1]},"2"={PCR[,1]=PC_INV2[,1]},"3"={PCR[,1]=PCr[,2]},"4"={PCR[,1]=PC_INV2[,2]},"5"={PCR[,1]=PC[,1]},"6"={PCR[,1]=PC_INV1[,1]},"7"={PCR[,1]=PC[,2]},"8"={PCR[,1]=PC_INV1[,2]})
      switch(BB,"1"={PCR[,2]=PCr[,1]},"2"={PCR[,2]=PC_INV2[,1]},"3"={PCR[,2]=PCr[,2]},"4"={PCR[,2]=PC_INV2[,2]},"5"={PCR[,2]=PC[,1]},"6"={PCR[,2]=PC_INV1[,1]},"7"={PCR[,2]=PC[,2]},"8"={PCR[,2]=PC_INV1[,2]})
      
      }
      
      #OLS
      resid<-matrix(0,t,r)
      ar<-matrix(0,r,1)
      for (q in 1:r) {
        fit<-ar.ols(PCR[,q],aic=FALSE,order.max=1,demean=FALSE,intercept=FALSE)
        resid[,q]<-scale(fit$resid,center=TRUE,scale=TRUE)
        resid[is.na(resid)]<-mean(resid[-1,])
        ar[q]<-fit$ar
      }
  #Bootstrap AR
  for (bb in 1:B) {
    #Primer Bootstrap. Objetivo: obtener las Load Bootstrap
    PC_pseudo1<-matrix(0,t,r)
    resid_boot<-matrix(0,t,r)
    resamplet <- sample(nrow(Resid),size = N,replace = TRUE,prob = NULL)
    Resid_boot<-Resid[resamplet,]
    #Bootstrapeamos residuos del factor e idiosincrático para obtener nuevas Y
    for (q in 1:r) {
      resid_boot[,q]<-sample(resid[,q],size = t,replace = TRUE,prob = NULL)
      PC_pseudo1[1,q]=PCR[1,q]
      for (e in 2:t) {
        PC_pseudo1[e,q]<-ar[q]*PC_pseudo1[e-1,q]+sqrt(1-ar[q]^2)*resid_boot[e,q]                
      }
    }
     Yboot1<-Load%*%t(PC_pseudo1)+Resid_boot
    
    #BAI
    # Estimamos cargas y ar bootstrap
    Y<-t(Yboot1)
    RR<-Y%*%t(Y)
    eR<-eigen(RR)
    values<-eR$values[c(1:r)]
    vectors<-matrix(eR$vectors[,c(1:r)],t,r)
    PC_boot1<-sqrt(t)*vectors
    PC_boot1_inv<-(-1)*PC_boot1
    Load2<-matrix(t(Y)%*%PC_boot1/t,N,r)
    PC_boot1r<-PC_boot1%*%solve(solve(t(Load2)%*%Load2)%*%(t(Load2)%*%Load))
    PC_boot2_inv<-(-1)*PC_boot1r
    
    if (r==1) {
      aa<-lm(PC_pseudo1[,1] ~ 0+PC_boot1r[,1])
      bbb<-lm(PC_pseudo1[,1] ~ 0+PC_boot2_inv[,1])
      ee<-lm(PC_pseudo1[,1] ~ 0+PC_boot1[,1])
      ff<-lm(PC_pseudo1[,1] ~ 0+PC_boot1_inv[,1])
      
      
      
      
      RS1<-matrix(c(summary(aa)$coefficients[1, 1],summary(bbb)$coefficients[1, 1],summary(ee)$coefficients[1, 1],summary(ff)$coefficients[1, 1]),4,1)
                  
      AA = which.max(RS1)
      switch(AA,"1"={PC_boot1R[,1]=PC_boot1r[,1]},"2"={PC_boot1R[,1]=PC_boot2_inv[,1]},"3"={PC_boot1R[,1]=PC_boot1[,1]},"4"={PC_boot1R[,1]=PC_boot1_inv[,1]}) 
      
    } else {
    
    aa<-lm(PC_pseudo1[,1] ~ 0+PC_boot1r[,1])
    bbb<-lm(PC_pseudo1[,1] ~ 0+PC_boot2_inv[,1])
    cc<-lm(PC_pseudo1[,1] ~ 0+PC_boot1r[,2])
    dd<-lm(PC_pseudo1[,1] ~ 0+PC_boot2_inv[,2])
    ee<-lm(PC_pseudo1[,1] ~ 0+PC_boot1[,1])
    ff<-lm(PC_pseudo1[,1] ~ 0+PC_boot1_inv[,1])
    gg<-lm(PC_pseudo1[,1] ~ 0+PC_boot1[,2])
    hh<-lm(PC_pseudo1[,1] ~ 0+PC_boot1_inv[,2])
    
    ii<-lm(PC_pseudo1[,2] ~ 0+PC_boot1r[,1])
    jjj<-lm(PC_pseudo1[,2] ~ 0+PC_boot2_inv[,1])
    kk<-lm(PC_pseudo1[,2] ~ 0+PC_boot1r[,2])
    ll<-lm(PC_pseudo1[,2] ~ 0+PC_boot2_inv[,2])
    mm<-lm(PC_pseudo1[,2] ~ 0+PC_boot1[,1])
    nn<-lm(PC_pseudo1[,2] ~ 0+PC_boot1_inv[,1])
    oo<-lm(PC_pseudo1[,2] ~ 0+PC_boot1[,2])
    pp<-lm(PC_pseudo1[,2] ~ 0+PC_boot1_inv[,2])
    
    RS1<-matrix(c(summary(aa)$coefficients[1, 1],summary(bbb)$coefficients[1, 1],summary(cc)$coefficients[1, 1],summary(dd)$coefficients[1, 1],summary(ee)$coefficients[1, 1],summary(ff)$coefficients[1, 1],summary(gg)$coefficients[1, 1],summary(hh)$coefficients[1, 1]),8,1)
    RS2<-matrix(c(summary(ii)$coefficients[1, 1],summary(jjj)$coefficients[1, 1],summary(kk)$coefficients[1, 1],summary(ll)$coefficients[1, 1],summary(mm)$coefficients[1, 1],summary(nn)$coefficients[1, 1],summary(oo)$coefficients[1, 1],summary(pp)$coefficients[1, 1]),8,1)
    
    AA = which.max(RS1)
    BB= which.max(RS2)
    switch(AA,"1"={PC_boot1R[,1]=PC_boot1r[,1]},"2"={PC_boot1R[,1]=PC_boot2_inv[,1]},"3"={PC_boot1R[,1]=PC_boot1r[,2]},"4"={PC_boot1R[,1]=PC_boot2_inv[,2]},"5"={PC_boot1R[,1]=PC_boot1[,1]},"6"={PC_boot1R[,1]=PC_boot1_inv[,1]},"7"={PC_boot1R[,1]=PC_boot1[,2]},"8"={PC_boot1R[,1]=PC_boot1_inv[,2]})
    switch(BB,"1"={PC_boot1R[,2]=PC_boot1r[,1]},"2"={PC_boot1R[,2]=PC_boot2_inv[,1]},"3"={PC_boot1R[,2]=PC_boot1r[,2]},"4"={PC_boot1R[,2]=PC_boot2_inv[,2]},"5"={PC_boot1R[,2]=PC_boot1[,1]},"6"={PC_boot1R[,2]=PC_boot1_inv[,1]},"7"={PC_boot1R[,2]=PC_boot1[,2]},"8"={PC_boot1R[,2]=PC_boot1_inv[,2]})
    }
    Loadb<-matrix(t(Y)%*%PC_boot1R/t,N,r)
    ar_boot1<-matrix(0,r,1)
    
   
   
    
    
    for (q in 1:r) {
      fit2<-ar.ols(PC_boot1R[,q],aic=FALSE,order.max=1,demean=FALSE,intercept=FALSE)
      ar_boot1[q]<-fit2$ar
    }
    resamplet <- sample(nrow(Resid),size = N,replace = TRUE,prob = NULL)
    Resid_boot<-Resid[resamplet,]
    for (tt in 2:t) {
      PC_pseudo2[tt,]<-ar_boot1[q]*PCR[tt-1,]+sqrt(1-ar_boot1[q]^2)*resid[tt]
    }
    Y_boot2<-matrix(Loadb%*%t(PC_pseudo2)+Resid_boot,nrow=N,ncol=t)
  
#Sacando PC del bootstrap
    #BAI
    Y<-t(Y_boot2)
    RR<-Y%*%t(Y)
    eR<-eigen(RR)
    values<-eR$values[c(1:r)]
    vectors<-matrix(eR$vectors[,c(1:r)],t,r)
    PC_boot2<-sqrt(t)*vectors
   PC_boot2_inv<-(-1)*PC_boot2
   Load3<-matrix(t(Y)%*%PC_boot2/t,N,r)
    PC_boot2r<- PC_boot2#%*%solve(solve(t(Load)%*%Load)%*%(t(Load)%*%P))
   PC_boot2r_inv<-(-1)*PC_boot2r
if (r==1) {
  aa<-lm(PC_pseudo2[,1] ~ 0+PC_boot2r[,1])
  bbb<-lm(PC_pseudo2[,1] ~ 0+PC_boot2r_inv[,1])
  ee<-lm(PC_pseudo2[,1] ~ 0+PC_boot2[,1])
  ff<-lm(PC_pseudo2[,1] ~ 0+PC_boot2_inv[,1])
  
  
  
  
  RS1<-matrix(c(summary(aa)$coefficients[1, 1],summary(bbb)$coefficients[1, 1],summary(ee)$coefficients[1, 1],summary(ff)$coefficients[1, 1]),4,1)
              
  AA = which.max(RS1)
  switch(AA,"1"={PC_boot2R[,1]=PC_boot2r[,1]},"2"={PC_boot2R[,1]=PC_boot2r_inv[,1]},"3"={PC_boot2R[,1]=PC_boot2[,1]},"4"={PC_boot1R[,1]=PC_boot2_inv[,1]}) 
  
  PC1_BOOT[,bb]<-PC_boot2R[,1]             
} else {

aa<-lm(PC_pseudo2[,1] ~ 0+PC_boot2r[,1])
bbb<-lm(PC_pseudo2[,1] ~ 0+PC_boot2r_inv[,1])
cc<-lm(PC_pseudo2[,1] ~ 0+PC_boot2r[,2])
dd<-lm(PC_pseudo2[,1] ~ 0+PC_boot2r_inv[,2])
ee<-lm(PC_pseudo2[,1] ~ 0+PC_boot2[,1])
ff<-lm(PC_pseudo2[,1] ~ 0+PC_boot2_inv[,1])
gg<-lm(PC_pseudo2[,1] ~ 0+PC_boot2[,2])
hh<-lm(PC_pseudo2[,1] ~ 0+PC_boot2_inv[,2])

ii<-lm(PC_pseudo2[,2] ~ 0+PC_boot2r[,1])
jjj<-lm(PC_pseudo2[,2] ~ 0+PC_boot2r_inv[,1])
kk<-lm(PC_pseudo2[,2] ~ 0+PC_boot2r[,2])
ll<-lm(PC_pseudo2[,2] ~ 0+PC_boot2r_inv[,2])
mm<-lm(PC_pseudo2[,2] ~ 0+PC_boot2[,1])
nn<-lm(PC_pseudo2[,2] ~ 0+PC_boot2_inv[,1])
oo<-lm(PC_pseudo2[,2] ~ 0+PC_boot2[,2])
pp<-lm(PC_pseudo2[,2] ~ 0+PC_boot2_inv[,2])

RS1<-matrix(c(summary(aa)$coefficients[1, 1],summary(bbb)$coefficients[1, 1],summary(cc)$coefficients[1, 1],summary(dd)$coefficients[1, 1],summary(ee)$coefficients[1, 1],summary(ff)$coefficients[1, 1],summary(gg)$coefficients[1, 1],summary(hh)$coefficients[1, 1]),8,1)
RS2<-matrix(c(summary(ii)$coefficients[1, 1],summary(jjj)$coefficients[1, 1],summary(kk)$coefficients[1, 1],summary(ll)$coefficients[1, 1],summary(mm)$coefficients[1, 1],summary(nn)$coefficients[1, 1],summary(oo)$coefficients[1, 1],summary(pp)$coefficients[1, 1]),8,1)

AA = which.max(RS1)
BB= which.max(RS2)
switch(AA,"1"={PC_boot2R[,1]=PC_boot2r[,1]},"2"={PC_boot2R[,1]=PC_boot2r_inv[,1]},"3"={PC_boot2R[,1]=PC_boot2r[,2]},"4"={PC_boot2R[,1]=PC_boot2r_inv[,2]},"5"={PC_boot2R[,1]=PC_boot2[,1]},"6"={PC_boot2R[,1]=PC_boot2_inv[,1]},"7"={PC_boot2R[,1]=PC_boot2[,2]},"8"={PC_boot2R[,1]=PC_boot2_inv[,2]})
switch(BB,"1"={PC_boot2R[,2]=PC_boot2r[,1]},"2"={PC_boot2R[,2]=PC_boot2r_inv[,1]},"3"={PC_boot2R[,2]=PC_boot2r[,2]},"4"={PC_boot2R[,2]=PC_boot2r_inv[,2]},"5"={PC_boot2R[,2]=PC_boot2[,1]},"6"={PC_boot2R[,2]=PC_boot2_inv[,1]},"7"={PC_boot2R[,2]=PC_boot2[,2]},"8"={PC_boot2R[,2]=PC_boot2_inv[,2]})

PC1_BOOT[,bb]<-PC_boot2R[,1]
PC2_BOOT[,bb]<-PC_boot2R[,2]
}
par(mfrow=c(2,2))
ts.plot(F[1,])
lines(PCR[,1],col="blue") 
ts.plot(PC_pseudo1[,1])
lines(PC_boot1R[,1],col="red") 
ts.plot(F[1,])
lines(PCR[,1],col="blue") 
ts.plot(PC_pseudo2[,1])
lines(PC_boot2R[,1],col="red") 
  }
  
#Varianza Bootstrap, mse
MSE<-((t(F)-PCR)^2)
RMSEi<-matrix(sqrt(colMeans(MSE)),2,1)
var1_boot<-rowVars(PC1_BOOT)
desv1_boot<-sqrt(var1_boot)
Error1_boot<-desv1_boot
RMSEBI1<-sqrt(mean(var1_boot))

var2_boot<-rowVars(PC2_BOOT)
desv2_boot<-sqrt(var2_boot)
Error2_boot<-desv2_boot
RMSEBI2<-sqrt(mean(var2_boot))

RMSE1[jj]<-RMSEi[1,1]
RMSE2[jj]<-RMSEi[2,1]
RMSEB1[jj]<-RMSEBI1
RMSEB2[jj]<-RMSEBI2

#Confidence Intervals
for (qq in 1:t) {
#IC1
CIL1[qq]<-PCR[qq,1]-z2*Error1_boot[qq]
CIU1[qq]<-PCR[qq,1]+z2*Error1_boot[qq]
#IC2
#CIL2[qq]<-PCR[qq,2]-z2*Error2_boot[qq]
#CIU2[qq]<-PCR[qq,2]+z2*Error2_boot[qq]
}


#Coverages1
cover1[1] <- if (CIL1[1]<F[1,1]&CIU1[1]>F[1,1]) cover1[1]=0 else cover1[1]=1
for(zz in 2:t) {
  cover1[zz] <- if (CIL1[zz]<F[1,zz]&CIU1[zz]>F[1,zz]) cover1[zz]=0 else cover1[zz]=1
}
COVER1<-matrix(cover1,nrow=t,ncol=1)
COVERAGES1=COVERAGES1+COVER1

#coverages2
#cover2[1] <- if (CIL2[1]<F[1,1]&CIU2[1]>F[1,1]) cover2[1]=0 else cover2[1]=1
#for(zz in 2:t) {
 # cover2[zz] <- if (CIL2[zz]<F[1,zz]&CIU2[zz]>F[1,zz]) cover2[zz]=0 else cover2[zz]=1
#}
#COVER2<-matrix(cover2,nrow=t,ncol=1)
#COVERAGES2=COVERAGES2+COVER2



}
MRMSE1[rrr]<-colMeans(RMSE1)
MRMSE2[rrr]<-colMeans(RMSE2)
MRMSEB1[rrr]<-colMeans(RMSEB1)
MRMSEB2[rrr]<-colMeans(RMSEB2)
coverage_ratio1[rrr]=1-colSums(COVERAGES1)/(t*R)
#coverage_ratio2[rrr]=1-colSums(COVERAGES2)/(t*R)

www[rrr]<-SN
phi[rrr]<-phi1
}
}


salida <- data.frame(www,phi,MRMSE1,MRMSE2,MRMSEB1,MRMSEB2,coverage_ratio1, coverage_ratio2)
#write.csv2(salida,"M11 r=2 boot1 rotación propia B=300 T=100 N=200.csv")
