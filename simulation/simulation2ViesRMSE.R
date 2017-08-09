library(sn)
library(bssn)
library(ClusterR)
#install.packages("bssn")

setwd("C:/Users/luis_/Dropbox/Luis y Rocio/Research/Pacotes/bssn/R")

source("functions.bssmn.R")
source("functions.mixbssmn.R")
source("functions.R")
source("algEMmixbssmn.R")

pii                    <- c(0.6,0.4)
g                      <- 2

pii                    <- c(0.6,0.4)
g                      <- length(pii)
alpha                  <- c(0.25,0.25)
beta                   <- c(1,5)
mfamily                <- "t"
nu                     <- 3


replicate              <- 1000
n                      <- c(80,150,500,2500)

vies = rmse = rb       <- array(0,dim=c(length(n),ncol=6))
estimaciones           <- array(0,dim=c(replicate,ncol=6))
colnames(estimaciones) <- c("alpha1","alpha2","beta1","beta2","pii1","nu")
colnames(vies)         <- c("VIESalpha1","VIESalpha2","VIESbeta1","VIESbeta2","VIESpii1","VIESnu")
colnames(rmse)         <- c("RMSEalpha1","RMSElpha2","RMSEbeta1","RMSEbeta2","RMSEpii1","RMSEnu")
colnames(rb)           <- c("RBalpha1","RBalpha2","RBbeta1","RBbeta2","RBpii1","RBnu")
rownames(vies)         = rownames(rmse) = rownames(rb) <- c("80", "150","500", "2500")
vt                     <- c(alpha,beta, pii[1],nu)#valores teoricos
vtmatrix               <- matrix(c(vt),nc=ncol(rb),nr=replicate,byrow=TRUE)

start.timeI    <- proc.time()
for(k in 1:length(n))
{
  i                      <- 1
  while(i <= replicate)
  {
    ti                   <- rmixbssmn(n[k],alpha,beta,nu,pii, mfamily)$y; hist(ti,breaks = 25, freq = FALSE, main="");  sseq    <- seq(0.1,100,0.01); dens    <- d.mixed.bssmn(sseq, pii, alpha, beta, nu, mfamily);lines(sseq,dens,col="red",lwd=2)
    initial              <- try(initial.Values(ti,g=2,algorithm="k-means",mfamily,lower=0.1,upper=0.9,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))
    fit                  <- try(EMmixbssmn2(ti,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))
    if(class(fit)!="try-error" && class(initial)!="try-error" && fit$result$convergence == TRUE && fit$result$pii > 0.5 && fit$result$nu < 4)
    {
      print(i)
      estimaciones[i,]   <- c(fit$result$alpha,fit$result$beta,fit$result$pii[1],fit$result$nu)
      i                  <- i+1
    }
  }
  vies[k,]            <- (apply(estimaciones-vtmatrix,2,mean))
  rmse[k,]            <- sqrt(apply((estimaciones-vtmatrix)^2,2,mean))
  rb[k,]              <- abs(apply((estimaciones-vtmatrix)/vt,2,mean))
  #df                  <- data.frame(vies[k,],rmse[k,],rb[k,])
  #write.csv(round(df, digits=4), paste("FM_BS_VIES_RMSE_RB_resultadosSimul2","tamanho_amostra=",n[k],'replicate=',replicate,'csv',sep="."))
}
end.timeF <- proc.time() - start.timeI #Reset time
text      <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")
df        <- round(data.frame(vies,rmse,rb),digits = 4)
write.csv(df  , paste("FM-BS-resultadosSimul2VIES-RMSE-k-bumps Escenario 1_2",'csv',sep="."))



#Graphics
setwd("C:/Users/luis_/Dropbox/Luis y Rocio/Research/Pacotes/bssn2/Other/FM-BS-SMN/Simulaciones")
result1_simul2_n  <- read.csv(file="FM-BS-resultadosSimul2VIES-RMSE-k-bumps Escenario 1_2.csv",head=TRUE,sep=",")[,-1]
result1_simul2_t  <- read.csv(file="FM-BS-resultadosSimul2VIES-RMSE-k-bumps Escenario 1_2t.csv",head=TRUE,sep=",")
result1_simul2_sl <- read.csv(file="FM-BS-resultadosSimul2VIES-RMSE-k-bumps Escenario 1_2sl.csv",head=TRUE,sep=",")
result1_simul2_cn <- read.csv(file="FM-BS-resultadosSimul2VIES-RMSE-k-bumps Escenario 1_2cn.csv",head=TRUE,sep=",")

n=c(80,150,500,2500)

#Graficos  #dev.size("in") #comando para conhecer o tamanho da tela
gray1 <- rgb(51,51,51,max=255)
gray2 <- rgb(112,112,112,max=255)
gray3 <- rgb(150,150,150,max=255)


###########################################################################################
postscript("expbiasAlpha1joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(-0.03,0.004),type="n",xlab="Samples Sizes (n)", ylab="Bias", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(-0.03,0.004,0.01))
s1 <- smooth.spline(n,result1_simul2_n[,2] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,2] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,2] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,2] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(alpha)[1]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("bottomright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expbiasAlpha2joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(-0.03,0.004),type="n",xlab="Samples Sizes (n)", ylab="Bias", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(-0.03,0.004,0.01))
s1 <- smooth.spline(n,result1_simul2_n[,3] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,3] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,3] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,3] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(alpha)[2]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("bottomright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expbiasBeta1joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(-0.03,0.3),type="n",xlab="Samples Sizes (n)", ylab="Bias", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(-0.03,0.3,0.1))
s1 <- smooth.spline(n,result1_simul2_n[,4] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,4] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,4] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,4] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(beta)[1]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expbiasBeta2joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(-0.2,0.02),type="n",xlab="Samples Sizes (n)", ylab="Bias", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(-0.2,0.02,0.1))
s1 <- smooth.spline(n,result1_simul2_n[,5] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,5] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,5] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,5] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(beta)[2]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("bottomright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expbiasP1joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(-0.01,0.005),type="n",xlab="Samples Sizes (n)", ylab="Bias", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(-0.01,0.005,0.005))
s1 <- smooth.spline(n,result1_simul2_n[,6] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,6] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,6] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,6] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(p)[1]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("bottomright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript
###########################################################################################


###########################################################################################
postscript("expRMSEAlpha1joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(0,0.08),type="n",xlab="Samples Sizes (n)", ylab="RMSE", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(0,0.08,0.02))
s1 <- smooth.spline(n,result1_simul2_n[,7] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,8] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,8] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,9] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(alpha)[1]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expRMSEAlpha2joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(0,0.06),type="n",xlab="Samples Sizes (n)", ylab="RMSE", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(0,0.06,0.01))
s1 <- smooth.spline(n,result1_simul2_n[,8] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,9] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,9] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,10] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(alpha)[2]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expRMSEBeta1joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(0,0.9),type="n",xlab="Samples Sizes (n)", ylab="RMSE", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(0,0.9,0.1))
s1 <- smooth.spline(n,result1_simul2_n[,9] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,10] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,10] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,11] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(beta)[1]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expRMSEBeta2joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(0,1.6),type="n",xlab="Samples Sizes (n)", ylab="RMSE", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(0,1.6,0.4))
s1 <- smooth.spline(n,result1_simul2_n[,10] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,11] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,11] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,12] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(beta)[2]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("expRMSEP1joint.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(c(0,2500),c(0,0.06),type="n",xlab="Samples Sizes (n)", ylab="RMSE", xaxt="n", yaxt="n")
axis(1,seq(0,2500,500))
axis(2,seq(0,0.06,0.01))
s1 <- smooth.spline(n,result1_simul2_n[,11] , spar=0.25)
s2 <- smooth.spline(n,result1_simul2_t[,12] , spar=0.25)
s3 <- smooth.spline(n,result1_simul2_sl[,12] , spar=0.25)
s4 <- smooth.spline(n,result1_simul2_cn[,13] , spar=0.25)
lines(s1,type="o",pch=15,col="red",lwd=2,lty=2)
lines(s2,type="o",pch=16,col="blue",lwd=2,lty=2)
lines(s3,type="o",pch=16,col="black",lwd=2,lty=1)
lines(s4,type="o",pch=16,col="green",lwd=2,lty=2)
title(expression(paste(p)[1]))   #Titulo do grafico
abline(h=0,lty=3,lwd=1)   #Linha horizontal
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("red","blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript
###########################################################################################
