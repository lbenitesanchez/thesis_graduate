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
setwd("C:/Users/luis_/Dropbox/Luis y Rocio/Research/Pacotes/bssn/R")

source("functions.bssmn.R")
source("functions.mixbssmn.R")
source("functions.R")
source("algEMmixbssmn.R")
#source("~/Dropbox/Pacotes/bssn/R/Other/FM-BS-SMN/khmeans.R")
#source("C:/Users/Larissa/Dropbox/Pacotes/bssn/R/Other/FM-BS-SMN/khmeans.R")

n                      <- 1000
pii                    <- c(0.6,0.4)
g                      <- length(pii)
alpha                  <- c(0.25,0.25)
beta                   <- c(1,5)

replicate              <- 500

################################################################################################
#                                        Familia "t"
################################################################################################
mfamily                <- "t"
nu                     <- 3

estimaciones           <- matrix(0,nrow=replicate,ncol=3*g) #t and sl
EP                     <- matrix(0,nrow=replicate,ncol=3*g - 1)
colnames(estimaciones) <- c("alpha1","alpha2","beta1","beta2","pii1","nu")
colnames(EP)           <- c("EPalpha1","EPalpha2","EPbeta1","EPbeta2","EPpii1")

i <- 1
start.time             <- Sys.time()
while(i <= replicate)
{
  ti                   <- rmixbssmn(n,alpha,beta,nu,pii, mfamily)$y; hist(ti,breaks = 20, freq = FALSE, main="");  sseq    <- seq(0.1,20,0.01); dens    <- d.mixed.bssmn(sseq, pii, alpha, beta, nu, mfamily);lines(sseq,dens,col="red",lwd=2)
  initial              <- try(initial.Values(ti,g=2,algorithm="k-mean",mfamily,lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))
  fit                  <- try(algEMmixbssmn(ti,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))
  fit                  <- EMmixbssmn(ti,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE)
  EPt                  <- try(InfmatrixmixBSSMN(ti,fit$result$pii,fit$result$alpha,fit$result$beta,fit$result$nu,mfamily)$EP)
  if(class(fit)!="try-error" &&  class(initial)!="try-error" && class(EPt)!="try-error" && fit$result$convergence == TRUE && fit$result$pii[1] > 0.5 && fit$result$nu < 4)
  {
    cat(i,"to",replicate,"\n")
    estimaciones[i,]   <- c(fit$result$alpha,fit$result$beta,fit$result$pii[1],fit$result$nu)
    EP[i,]             <- EPt
    i                  <- i+1
  }
}
end.time               <- Sys.time()
time.taken             <- end.time - start.time

#+----------------------------------------------------------------------------------------------------+
#                                   Obtencao da COV FM-BS-SMN                                         #
#+----------------------------------------------------------------------------------------------------+

alpha1_csup            <- estimaciones[,1] + 1.96*EP[,1]; alpha1_cinf    <- estimaciones[,1] - 1.96*EP[,1]
alpha2_csup            <- estimaciones[,2] + 1.96*EP[,2]; alpha2_cinf    <- estimaciones[,2] - 1.96*EP[,2]
count.a1               <- c()
count.a2               <- c()

beta1_csup             <- estimaciones[,3] + 1.96*EP[,3]; beta1_cinf    <- estimaciones[,3] - 1.96*EP[,3]
beta2_csup             <- estimaciones[,4] + 1.96*EP[,4]; beta2_cinf    <- estimaciones[,4] - 1.96*EP[,4]
count.b1               <- c()
count.b2               <- c()

pii1_csup              <- estimaciones[,5] + 1.96*EP[,5];  pii1_cinf  <- estimaciones[,5] - 1.96*EP[,5]
count.pii1             <- c()

for(i in 1:replicate)
{
  if(alpha1_cinf[i]  <= alpha[1] && alpha1_csup[i] >= alpha[1])   count.a1[i] <- 1 else count.a1[i] <- 0
  if(alpha2_cinf[i]  <= alpha[2] && alpha2_csup[i] >= alpha[2])   count.a2[i] <- 1 else count.a2[i] <- 0

  if(beta1_cinf[i]  <= beta[1] && beta1_csup[i] >= beta[1])   count.b1[i] <- 1 else count.b1[i] <- 0
  if(beta2_cinf[i]  <= beta[2] && beta2_csup[i] >= beta[2])   count.b2[i] <- 1 else count.b2[i] <- 0

  if(pii1_cinf[i]   <= pii[1]   && pii1_csup[i] >= pii[1]) count.pii1[i] <- 1 else count.pii1[i] <- 0
}

COVMC_alpha1         <- sum(count.a1)/replicate
COVMC_alpha2         <- sum(count.a2)/replicate

COVMC_beta1          <- sum(count.b1)/replicate
COVMC_beta2          <- sum(count.b2)/replicate

COVMC_pii1           <- sum(count.pii1)/replicate
COVT                 <- 100*c(COVMC_alpha1,COVMC_alpha2,COVMC_beta1,COVMC_beta2,COVMC_pii1,0)

teoricos             <- matrix(c(alpha,beta,pii[1],nu),nrow=1)
colnames(teoricos)   <- c("alpha1","alpha2","beta1","beta2","pii1","nu")
resultados <- cbind(c(teoricos),c(round(apply(estimaciones,2,mean), digits= 4)),c(round(apply(EP,2,mean), digits= 4),0),c(round(apply(estimaciones,2,sd), digits= 4)),c(COVT))
colnames(resultados) <- c("Teoricos","Estimativas","EP","MC","COV")
resultados

################################################################################################
#                                   Familia Normal Contaminada
################################################################################################
mfamily                <- "cn"
nu                     <- c(0.1,0.1)

estimaciones           <- matrix(0,nrow=replicate,ncol=3*g + 1) #cn
EP                     <- matrix(0,nrow=replicate,ncol=3*g - 1)
colnames(estimaciones) <- c("alpha1","alpha2","beta1","beta2","pii1","nu","gama")
colnames(EP)           <- c("EPalpha1","EPalpha2","EPbeta1","EPbeta2","EPpii1")

i <- 0
start.time             <- Sys.time()
k <- 0
while(i <= replicate)
{
  ti                   <- rmixbssmn(n,alpha,beta,nu,pii, mfamily)$y; hist(ti,breaks = 25, freq = FALSE, main="");  sseq    <- seq(0.1,100,0.01); dens    <- d.mixed.bssmn(sseq, pii, alpha, beta, nu, mfamily);lines(sseq,dens,col="red",lwd=2)
  initial              <- try(initial.Values(ti,g=2,algorithm="k-means",mfamily,lower=0.1,upper=0.9,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))
  fit                  <- try(EMmixbssmn2(ti,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))
  k                    <- k + 1
  EPt                  <- try(InfmatrixmixBSSMN(ti,fit$result$pii,fit$result$alpha,fit$result$beta,fit$result$nu,mfamily)$EP)
  if(class(fit)!="try-error" && class(initial)!="try-error" && class(EPt)!="try-error" && fit$result$convergence == TRUE && fit$result$pii > 0.5 && fit$result$nu < 4)
  {
    i                  <- i+1
    print(i)
    estimaciones[i,]   <- c(fit$result$alpha,fit$result$beta,fit$result$pii[1],fit$result$nu)
    EP[i,]             <- EPt
  }
}
end.time               <- Sys.time()
time.taken             <- end.time - start.time


#+----------------------------------------------------------------------------------------------------+
#                                   Obtencao da COV FM-BS                                             #
#+----------------------------------------------------------------------------------------------------+

alpha1_csup            <- estimaciones[,1] + 1.96*EP[,1]; alpha1_cinf    <- estimaciones[,1] - 1.96*EP[,1]
alpha2_csup            <- estimaciones[,2] + 1.96*EP[,2]; alpha2_cinf    <- estimaciones[,2] - 1.96*EP[,2]
count.a1               <- c()
count.a2               <- c()

beta1_csup             <- estimaciones[,3] + 1.96*EP[,3]; beta1_cinf    <- estimaciones[,3] - 1.96*EP[,3]
beta2_csup             <- estimaciones[,4] + 1.96*EP[,4]; beta2_cinf    <- estimaciones[,4] - 1.96*EP[,4]
count.b1               <- c()
count.b2               <- c()

pii1_csup              <- estimaciones[,5] + 1.96*EP[,5];  pii1_cinf  <- estimaciones[,5] - 1.96*EP[,5]
count.pii1             <- c()

for(i in 1:replicate)
{
  if(alpha1_cinf[i]  <= alpha[1] && alpha1_csup[i] >= alpha[1])   count.a1[i] <- 1 else count.a1[i] <- 0
  if(alpha2_cinf[i]  <= alpha[2] && alpha2_csup[i] >= alpha[2])   count.a2[i] <- 1 else count.a2[i] <- 0

  if(beta1_cinf[i]  <= beta[1] && beta1_csup[i] >= beta[1])   count.b1[i] <- 1 else count.b1[i] <- 0
  if(beta2_cinf[i]  <= beta[2] && beta2_csup[i] >= beta[2])   count.b2[i] <- 1 else count.b2[i] <- 0

  if(pii1_cinf[i]   <= pii[1]   && pii1_csup[i] >= pii[1]) count.pii1[i] <- 1 else count.pii1[i] <- 0
}

COVMC_alpha1            <- sum(count.a1)/replicate
COVMC_alpha2            <- sum(count.a2)/replicate

COVMC_beta1             <- sum(count.b1)/replicate
COVMC_beta2             <- sum(count.b2)/replicate

COVMC_pii1              <- sum(count.pii1)/replicate
COVT                    <- 100*c(COVMC_alpha1,COVMC_alpha2,COVMC_beta1,COVMC_beta2,COVMC_pii1)



teoricos             <- matrix(c(alpha,beta,pii[1],nu),nrow=1)
colnames(teoricos)   <- c("alpha1","alpha2","beta1","beta2","pii1","nu","gamma")
resultados <- cbind(c(teoricos),c(round(apply(estimaciones,2,mean), digits= 4)),c(round(apply(EP,2,mean), digits= 4),0,0),c(round(apply(estimaciones,2,sd), digits= 4)),c(COVT,0,0))
colnames(resultados) <- c("Teoricos","Estimativas","EP","MC","COV")
resultados

################################################################################################
#                                        Familia "sl"
################################################################################################
mfamily                <- "sl"
nu                     <- 3

estimaciones           <- matrix(0,nrow=replicate,ncol=3*g) #t and sl
EP                     <- matrix(0,nrow=replicate,ncol=3*g - 1)
colnames(estimaciones) <- c("alpha1","alpha2","beta1","beta2","pii1","nu")
colnames(EP)           <- c("EPalpha1","EPalpha2","EPbeta1","EPbeta2","EPpii1")

i <- 1
start.time             <- Sys.time()
while(i <= replicate)
{
  ti                   <- rmixbssmn(n,alpha,beta,nu,pii, mfamily)$y; hist(ti,breaks = 20, freq = FALSE, main="");  sseq    <- seq(0.1,20,0.01); dens    <- d.mixed.bssmn(sseq, pii, alpha, beta, nu, mfamily);lines(sseq,dens,col="red",lwd=2)
  initial              <- try(initial.Values(ti,g=2,algorithm="k-medoids",mfamily,lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE))
  fit                  <- try(EMmixbssmn2(ti,initial$alpha, initial$beta, initial$nu, initial$pii, g , mfamily, accuracy = 10^-6, iter.max = 500,aitken=TRUE))
  EPt                  <- try(InfmatrixmixBSSMN(ti,fit$result$pii,fit$result$alpha,fit$result$beta,fit$result$nu,mfamily)$EP)
  if(class(fit)!="try-error" &&  class(initial)!="try-error" && class(EPt)!="try-error" && fit$result$convergence == TRUE && fit$result$pii[1] > 0.5 && fit$result$nu < 4)
  {
    cat(i,"to",replicate,"\n")
    estimaciones[i,]   <- c(fit$result$alpha,fit$result$beta,fit$result$pii[1],fit$result$nu)
    EP[i,]             <- EPt
    i                  <- i+1
  }
}
end.time               <- Sys.time()
time.taken             <- end.time - start.time

#+----------------------------------------------------------------------------------------------------+
#                                   Obtencao da COV FM-BS-SMN                                         #
#+----------------------------------------------------------------------------------------------------+

alpha1_csup            <- estimaciones[,1] + 1.96*EP[,1]; alpha1_cinf    <- estimaciones[,1] - 1.96*EP[,1]
alpha2_csup            <- estimaciones[,2] + 1.96*EP[,2]; alpha2_cinf    <- estimaciones[,2] - 1.96*EP[,2]
count.a1               <- c()
count.a2               <- c()

beta1_csup             <- estimaciones[,3] + 1.96*EP[,3]; beta1_cinf    <- estimaciones[,3] - 1.96*EP[,3]
beta2_csup             <- estimaciones[,4] + 1.96*EP[,4]; beta2_cinf    <- estimaciones[,4] - 1.96*EP[,4]
count.b1               <- c()
count.b2               <- c()

pii1_csup              <- estimaciones[,5] + 1.96*EP[,5];  pii1_cinf  <- estimaciones[,5] - 1.96*EP[,5]
count.pii1             <- c()

for(i in 1:replicate)
{
  if(alpha1_cinf[i]  <= alpha[1] && alpha1_csup[i] >= alpha[1])   count.a1[i] <- 1 else count.a1[i] <- 0
  if(alpha2_cinf[i]  <= alpha[2] && alpha2_csup[i] >= alpha[2])   count.a2[i] <- 1 else count.a2[i] <- 0

  if(beta1_cinf[i]  <= beta[1] && beta1_csup[i] >= beta[1])   count.b1[i] <- 1 else count.b1[i] <- 0
  if(beta2_cinf[i]  <= beta[2] && beta2_csup[i] >= beta[2])   count.b2[i] <- 1 else count.b2[i] <- 0

  if(pii1_cinf[i]   <= pii[1]   && pii1_csup[i] >= pii[1]) count.pii1[i] <- 1 else count.pii1[i] <- 0
}

COVMC_alpha1         <- sum(count.a1)/replicate
COVMC_alpha2         <- sum(count.a2)/replicate

COVMC_beta1          <- sum(count.b1)/replicate
COVMC_beta2          <- sum(count.b2)/replicate

COVMC_pii1           <- sum(count.pii1)/replicate
COVT                 <- 100*c(COVMC_alpha1,COVMC_alpha2,COVMC_beta1,COVMC_beta2,COVMC_pii1,0)

teoricos             <- matrix(c(alpha,beta,pii[1],nu),nrow=1)
colnames(teoricos)   <- c("alpha1","alpha2","beta1","beta2","pii1","nu")
resultados <- cbind(c(teoricos),c(round(apply(estimaciones,2,mean), digits= 4)),c(round(apply(EP,2,mean), digits= 4),0),c(round(apply(estimaciones,2,sd), digits= 4)),c(COVT))
colnames(resultados) <- c("Teoricos","Estimativas","EP","MC","COV")
resultados



###################################
postscript("boxplotAlpha1Simul2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
boxplot(estimaciones[,1],main=expression(paste(alpha[1],"=0.25")))
dev.off() #Fechando o dispositivo potscript

postscript("boxplotAlpha2Simul2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
boxplot(estimaciones[,2],main=expression(paste(alpha[2],"=0.25")))
dev.off() #Fechando o dispositivo potscript

postscript("boxplotBeta1Simul2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
boxplot(estimaciones[,3],main=expression(paste(beta[1],"=1")))
dev.off() #Fechando o dispositivo potscript

postscript("boxplotBeta2Simul2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
boxplot(estimaciones[,4],main=expression(paste(beta[2],"=5")))
dev.off() #Fechando o dispositivo potscript
