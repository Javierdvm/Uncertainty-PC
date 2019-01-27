library('depth')
library('geometry')
####### Depth
depthtukey = vector(mode = "numeric", length = 0)
sco<-matrix(rnorm(5000,0,1),2500,2)
for (k in c(1:dim(sco)[1])) {
  depthtukey[k] = depth(sco[k,], sco, method ='Tukey')
}
tukey <- sort(depthtukey, decreasing = TRUE, index.return = TRUE); Jordered <- sco[tukey$ix]

plot(sco)
points(sco[tukey$ix[1],1] , sco[tukey$ix[1],2], col='red', lwd=3)

factor=0.95
points(sco[tukey$ix[1:length(tukey$x)*factor],1] , sco[tukey$ix[1:length(tukey$x)*factor],2] , col='red', lwd=3)

ScoFactor<-sco[tukey$ix[1:length(tukey$x)*factor],]
points<-t(convhulln(ScoFactor))
fuckingpoints<-cbind(ScoFactor[points,1],ScoFactor[points,2])

points(fuckingpoints[,1],fuckingpoints[,2], col='blue', lwd=3)

xnew <- fuckingpoints[,1][order(Arg(scale(fuckingpoints[,1]) + scale(fuckingpoints[,2]) * 1i))]
ynew <- fuckingpoints[,2][order(Arg(scale(fuckingpoints[,1]) + scale(fuckingpoints[,2]) * 1i))]

polygon(xnew ,ynew)
