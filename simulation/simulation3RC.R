# +---------------------------------------------------------------------------------------+ #
# +--------------------------    Mixture BS-SMN     --------------------------------------+ #
# +---------------------------------------------------------------------------------------+ #

library(sn)
library(bssn)
library(mixsmsn)
library(ClusterR) #Este pacote e para poder utilizar a algoritmo k-medoids

#Windows
#setwd("C:/Users/Larissa/Dropbox/Pacotes/bssn/R")
#setwd("C:/Users/Luiz/Dropbox/Pacotes/bssn/R")
setwd("C:/Users/Luiz/Dropbox/Luis y Rocio/Research/Pacotes/bssn/R")
setwd("C:/Users/luis_/Dropbox/Luis y Rocio/Research/Pacotes/bssn/R")

source("functions.bssmn.R")
source("functions.mixbssmn.R")
source("functions.R")
source("algEMmixbssmn.R")
source("functions.mixbs.R")
source("algEMmixbs.R")
#source("~/Dropbox/Pacotes/bssn/R/Other/FM-BS-SMN/khmeans.R")
#source("C:/Users/Larissa/Dropbox/Pacotes/bssn/R/Other/FM-BS-SMN/khmeans.R")

n                      <- 1000
pii                    <- c(0.6,0.4)
g                      <- length(pii)
alpha                  <- c(0.25,0.25)
beta                   <- c(1,5)

replicate              <- 1000
perturb = seq(1:10); obsp=150

################################################################################################
#                                        Caso FM-BS
################################################################################################

estimBSSMSN           <- array(0,dim=c(replicate,ncol=6))
estimBSSMSNp          <- array(0,dim=c(replicate,ncol=6))
RC                    <- array(0,dim=c(length(perturb),ncol=6))
colnames(RC)          <- c("RCAlpha1","RCAlpha2","RCBeta1","RCBeta2","RCpii1","RCpii2")

start.timeI           <- proc.time();
for(kk in 1:length(perturb))
{


  i                    <- 1
  while(i <= replicate)
  {#Begin While replicate
    #x                    <- gen.BS.Skew.t(n, delta, nu, alpha, beta1)
    #x                    <-gen.BS.Skew.slash(n, delta, nu, alpha, beta1)
    x                    <- rmixbs(n,alpha,beta,pii)$y;hist(x,breaks=30)
    xp                   <- x
    xp[obsp]             <- x[obsp] + perturb[kk]

    cat("#Replicate",i,"from a total of",replicate,",Simulation Time -->",round((proc.time() - start.timeI)[3]/60, digits=3),"minutes",'\n')

    est <- try(alg.EM.mix.bs(x, alpha, beta, pii, g=2, get.init = FALSE, accuracy = 10^-6, iter.max = 500, kmeans.param = NULL,aitken = TRUE))

    estp<- try(alg.EM.mix.bs(xp,  alpha, beta, pii, g=2, get.init = FALSE, accuracy = 10^-6, iter.max = 500, kmeans.param = NULL,aitken = TRUE))

    if(class(est)!="try-error" && est$result$convergence==TRUE && class(estp)!="try-error" && estp$result$convergence==TRUE)
    {
      estimBSSMSN[i,]    <- c(est$result$alpha,est$result$beta,est$result$pii)
      estimBSSMSNp[i,]   <- c(estp$result$alpha,estp$result$beta,estp$result$pii)
      i                  <- i+1
    }
  }#End While replicate
  RC[kk,]               <- abs((apply(estimBSSMSNp,2,mean) - apply(estimBSSMSN,2,mean))/apply(estimBSSMSN,2,mean))
}

#Save results
write.csv(RC  , paste("resultadosSimul3-mixbs",'csv',sep="."))

end.timeF             <- proc.time() - start.timeI #Reset time
text                  <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")
return(text)



################################################################################################
#                                        Familia "t"
################################################################################################
mfamily                <- "t"
nu                     <- 3

estimBSSMSN           <- array(0,dim=c(replicate,ncol=6))
estimBSSMSNp          <- array(0,dim=c(replicate,ncol=6))
RC                    <- array(0,dim=c(length(perturb),ncol=6))
colnames(RC)          <- c("RCAlpha1","RCAlpha2","RCBeta1","RCBeta2","RCpii1","RCpii2")

start.timeI           <- proc.time();
for(kk in 1:length(perturb))
{


  i                    <- 1
  while(i <= replicate)
  {#Begin While replicate
    #x                    <- gen.BS.Skew.t(n, delta, nu, alpha, beta1)
    #x                    <-gen.BS.Skew.slash(n, delta, nu, alpha, beta1)
    x                    <- rmixbssmn(n,alpha,beta,nu,pii, mfamily)$y; hist(x,breaks = 20, freq = FALSE, main="");  sseq    <- seq(0.1,20,0.01); dens    <- d.mixed.bssmn(sseq, pii, alpha, beta, nu, mfamily);lines(sseq,dens,col="red",lwd=2)
    xp                   <- x
    xp[obsp]             <- x[obsp] + perturb[kk]

    cat("#Replicate",i,"from a total of",replicate,",Simulation Time -->",round((proc.time() - start.timeI)[3]/60, digits=3),"minutes",'\n')

    initial              <- try(initial.Values(x,g=2,algorithm="k-means",mfamily,lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))
    initial1              <- try(initial.Values(xp,g=2,algorithm="k-means",mfamily,lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))

     est <- try(EMmixbssmn2(x,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))

    estp<- try(EMmixbssmn2(xp,initial1$alpha, initial1$beta, initial1$nu, initial1$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))

    if(class(est)!="try-error" && est$result$convergence==TRUE && class(estp)!="try-error" && estp$result$convergence==TRUE)
    {
      estimBSSMSN[i,]    <- c(est$result$alpha,est$result$beta,est$result$pii)
      estimBSSMSNp[i,]   <- c(estp$result$alpha,estp$result$beta,estp$result$pii)
      i                  <- i+1
    }
  }#End While replicate
  RC[kk,]               <- abs((apply(estimBSSMSNp,2,mean) - apply(estimBSSMSN,2,mean))/apply(estimBSSMSN,2,mean))
}

#Save results
write.csv(RC  , paste("resultadosSimul3t-mixbssmn",'csv',sep="."))

end.timeF             <- proc.time() - start.timeI #Reset time
text                  <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")
return(text)

################################################################################################
#                                   Familia Normal Contaminada
################################################################################################
mfamily                <- "cn"
nu                     <- c(0.1,0.1)

estimBSSMSN           <- array(0,dim=c(replicate,ncol=6))
estimBSSMSNp          <- array(0,dim=c(replicate,ncol=6))
RC                    <- array(0,dim=c(length(perturb),ncol=6))
colnames(RC)          <- c("RCAlpha1","RCAlpha2","RCBeta1","RCBeta2","RCpii1","RCpii2")

start.timeI           <- proc.time();
for(kk in 1:length(perturb))
{


  i                    <- 1
  while(i <= replicate)
  {#Begin While replicate
    #x                    <- gen.BS.Skew.t(n, delta, nu, alpha, beta1)
    #x                    <-gen.BS.Skew.slash(n, delta, nu, alpha, beta1)
    x                    <- rmixbssmn(n,alpha,beta,nu,pii, mfamily)$y; hist(x,breaks = 20, freq = FALSE, main="");  sseq    <- seq(0.1,20,0.01); dens    <- d.mixed.bssmn(sseq, pii, alpha, beta, nu, mfamily);lines(sseq,dens,col="red",lwd=2)
    xp                   <- x
    xp[obsp]             <- x[obsp] + perturb[kk]

    cat("#Replicate",i,"from a total of",replicate,",Simulation Time -->",round((proc.time() - start.timeI)[3]/60, digits=3),"minutes",'\n')

    initial              <- try(initial.Values(x,g=2,algorithm="k-means",mfamily,lower=0.1,upper=0.9,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))
    initial1              <- try(initial.Values(xp,g=2,algorithm="k-means",mfamily,lower=0.1,upper=0.9,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))

    est <- try(EMmixbssmn2(x,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))

    estp<- try(EMmixbssmn2(xp,initial1$alpha, initial1$beta, initial1$nu, initial1$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))

    if(class(est)!="try-error" && est$result$convergence==TRUE && class(estp)!="try-error" && estp$result$convergence==TRUE)
    {
      estimBSSMSN[i,]    <- c(est$result$alpha,est$result$beta,est$result$pii)
      estimBSSMSNp[i,]   <- c(estp$result$alpha,estp$result$beta,estp$result$pii)
      i                  <- i+1
    }
  }#End While replicate
  RC[kk,]               <- abs((apply(estimBSSMSNp,2,mean) - apply(estimBSSMSN,2,mean))/apply(estimBSSMSN,2,mean))
}

#Save results
write.csv(RC  , paste("resultadosSimul3cn-mixbssmn",'csv',sep="."))

end.timeF             <- proc.time() - start.timeI #Reset time
text                  <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")
return(text)

################################################################################################
#                                        Familia "sl"
################################################################################################
mfamily                <- "sl"
nu                     <- 3

estimBSSMSN           <- array(0,dim=c(replicate,ncol=6))
estimBSSMSNp          <- array(0,dim=c(replicate,ncol=6))
RC                    <- array(0,dim=c(length(perturb),ncol=6))
colnames(RC)          <- c("RCAlpha1","RCAlpha2","RCBeta1","RCBeta2","RCpii1","RCpii2")

start.timeI           <- proc.time();
for(kk in 1:length(perturb))
{


  i                    <- 1
  while(i <= replicate)
  {#Begin While replicate
    #x                    <- gen.BS.Skew.t(n, delta, nu, alpha, beta1)
    #x                    <-gen.BS.Skew.slash(n, delta, nu, alpha, beta1)
    x                    <- rmixbssmn(n,alpha,beta,nu,pii, mfamily)$y; hist(x,breaks = 20, freq = FALSE, main="");  sseq    <- seq(0.1,20,0.01); dens    <- d.mixed.bssmn(sseq, pii, alpha, beta, nu, mfamily);lines(sseq,dens,col="red",lwd=2)
    xp                   <- x
    xp[obsp]             <- x[obsp] + perturb[kk]

    cat("#Replicate",i,"from a total of",replicate,",Simulation Time -->",round((proc.time() - start.timeI)[3]/60, digits=3),"minutes",'\n')

    initial              <- try(initial.Values(x,g=2,algorithm="k-means",mfamily,lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))
    initial1              <- try(initial.Values(xp,g=2,algorithm="k-means",mfamily,lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))

    est <- try(EMmixbssmn2(x,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))

    estp<- try(EMmixbssmn2(xp,initial1$alpha, initial1$beta, initial1$nu, initial1$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))

    if(class(est)!="try-error" && est$result$convergence==TRUE && class(estp)!="try-error" && estp$result$convergence==TRUE)
    {
      estimBSSMSN[i,]    <- c(est$result$alpha,est$result$beta,est$result$pii)
      estimBSSMSNp[i,]   <- c(estp$result$alpha,estp$result$beta,estp$result$pii)
      i                  <- i+1
    }
  }#End While replicate
  RC[kk,]               <- abs((apply(estimBSSMSNp,2,mean) - apply(estimBSSMSN,2,mean))/apply(estimBSSMSN,2,mean))
}

#Save results
write.csv(RC  , paste("resultadosSimul3sl-mixbssmn",'csv',sep="."))

end.timeF             <- proc.time() - start.timeI #Reset time
text                  <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")
return(text)


setwd("C:/Users/luis_/Dropbox/Luis y Rocio/Research/Pacotes/bssn2/Other/FM-BS-SMN/Simulaciones")

result1=read.csv("resultadosSimul3normal-mixbssmn.csv",header=T)
result2=read.csv("resultadosSimul3t-mixbssmn.csv",header=T)
result3=read.csv("resultadosSimul3cn-mixbssmn.csv",header=T)
result4=read.csv("resultadosSimul3sl-mixbssmn.csv",header=T)

v                  <- c(1:10)
palette(gray(seq(0,.9,len =30)))
RCAlpha1           <- cbind(result4[,2],result3[,2],result2[,2],result1[,2])
RCAlpha2           <- cbind(result1[,3],result2[,3],result3[,3],result4[,3])
RCBeta1            <- cbind(result1[,4],result2[,4],result3[,4],result4[,4])
RCBeta2            <- cbind(result1[,5],result2[,5],result3[,5],result4[,5])
RCpii1             <- cbind(result1[,6],result2[,6],result3[,6],result4[,6])


colnames(RCAlpha1) <- c("MF-BS-N","MF-BS-T","MF-BS-NC","MF-BS-SL"); rownames(RCAlpha1) <- c(1:10)
colnames(RCAlpha2) <- c("MF-BS-N","MF-BS-T","MF-BS-NC","MF-BS-SL"); rownames(RCAlpha2) <- c(1:10)
colnames(RCBeta1)  <- c("MF-BS-N","MF-BS-T","MF-BS-NC","MF-BS-SL"); rownames(RCBeta1)  <- c(1:10)
colnames(RCBeta2)  <- c("MF-BS-N","MF-BS-T","MF-BS-NC","MF-BS-SL"); rownames(RCBeta2)  <- c(1:10)
colnames(RCpii1)   <- c("MF-BS-N","MF-BS-T","MF-BS-NC","MF-BS-SL"); rownames(RCpii1)   <- c(1:10)


postscript("exp3RCalpha1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(RCAlpha1, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.015), space=c(0,2))
abline(h=0.002, lty = 3, col="black")
abline(h=0.004, lty = 3, col="black")
abline(h=0.006, lty = 3, col="black")
abline(h=0.008, lty = 3, col="black")
abline(h=0.010, lty = 3, col="black")
abline(h=0.012, lty = 3, col="black")
box()
ypos1 <-  c(max(RCAlpha1)+0.0005, max(RCAlpha1)-0.008,max(RCAlpha1)-0.010,max(RCAlpha1)-0.0095)
text(x=aux[6,]-0.5,y=ypos1, label=expression(1---paste(vartheta)---10), cex=0.85)
title(expression(paste(alpha[1])))
dev.off() #Fechando o dispositivo potscript

postscript("exp3RCalpha2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(RCAlpha2, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.03), space=c(0,2))
abline(h=0.005, lty = 3, col="black")
abline(h=0.010, lty = 3, col="black")
abline(h=0.015, lty = 3, col="black")
box()
ypos1 <-  c(max(RCAlpha2)+0.002, max(RCAlpha2)-0.011,max(RCAlpha2)-0.013,max(RCAlpha2)-0.012)
text(x=aux[6,]-0.5,y=ypos1, label=expression(1---paste(vartheta)---10), cex=0.85)
title(expression(paste(alpha[2])))
dev.off() #Fechando o dispositivo potscript

postscript("exp3RCbeta1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(RCBeta1, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.001), space=c(0,2))
abline(h=0.0001, lty = 3, col="black")
abline(h=0.0004, lty = 3, col="black")
abline(h=0.0007, lty = 3, col="black")
box()
ypos1 <-  c(max(RCBeta1)+0.00005, max(RCBeta1)-0.0002, max(RCBeta1)-0.00025, max(RCBeta1)-0.0002)
text(x=aux[6,]-0.5,y=ypos1, label=expression(1---paste(vartheta)---10), cex=0.85)
title(expression(paste(beta[1])))
dev.off() #Fechando o dispositivo potscript

postscript("exp3RCbeta2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(RCBeta2, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.004), space=c(0,2))
abline(h=0.001, lty = 3, col="black")
abline(h=0.002, lty = 3, col="black")
abline(h=0.003, lty = 3, col="black")
box()
ypos1 <-  c(max(RCBeta2)+0.0002, max(RCBeta2)-0.0013, max(RCBeta2)-0.0011, max(RCBeta2)-0.0009)
text(x=aux[6,]-0.5,y=ypos1, label=expression(1---paste(vartheta)---10), cex=0.85)
title(expression(paste(beta[2])))
dev.off() #Fechando o dispositivo potscript

postscript("exp3RCpii1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(RCpii1, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.002), space=c(0,2))
abline(h=0.0005, lty = 3, col="black")
abline(h=0.0010, lty = 3, col="black")
abline(h=0.0015, lty = 3, col="black")
box()
ypos1 <-  c(max(RCpii1), max(RCpii1) + 0.0001, max(RCpii1) - 0.0001 , max(RCpii1))
text(x=aux[6,]-0.5,y=ypos1, label=expression(1---paste(vartheta)---10), cex=0.85)
title(expression(paste(p[1])))
dev.off() #Fechando o dispositivo potscript






