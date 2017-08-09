setwd("C:/Users/luis_/Dropbox/Luis y Rocio/Research/Pacotes/bssn/R")

source("functions.R")
source("functions.mixbssmn.R")

# Hazerd with 2 components
#First component
alpha1                 <- 1
beta1                  <- 1

#Second component
alpha2                 <- 2
beta2                  <- 2

alpha                  <- c(alpha1,alpha2)
beta                   <- c(beta1,beta2)
pii                    <- c(0.2,0.8)


postscript("hazarFM-BS1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
 c      <- 1/(2*alpha2^2*beta2)
 ti     <- seq(0.001,30,0.01)
 HMIXBS <- H.mix.bssmn(ti,pii,alpha,beta,mfamily="t")
 plot(ti,HMIXBS,type = "l",xlab="t", ylab="h(t)",ylim=c(0,1))
 abline(h=c,col="black",lty=2)
dev.off() #Fechando o dispositivo potscript

# ---------------------------------------------------------------------------------
# Hazerd with 1 component
#First component
alpha1                 <- 1
beta1                  <- 1

alpha                  <- c(alpha1)
beta                   <- c(beta1)
pii                    <- c(1)
#(alpha1/alpha2)^2*(beta1/beta2)<1#Condicion


postscript("hazarFM-BS1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
 c      <- 1/(2*alpha1^2*beta1)
 ti     <- seq(0.001,10,0.01)
 HMIXBS <- H.mix.bs(ti,pii,alpha,beta)
 plot(ti,HMIXBS,type = "l",xlab="t", ylab="h(t)",ylim=c(0,1))
 abline(h=c,col="red",lty=2)
dev.off() #Fechando o dispositivo potscript



#########################################################################
postscript("hazardgdif1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,40,0.01)
f1  <- H.mix.bs(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7))
f2  <- H.mix.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7),nu = 4, mfamily = "t")
f3  <- H.mix.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7),nu = 4, mfamily = "sl")
f4  <- H.mix.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "h(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("bottomleft",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("hazardgdif2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,25,0.01)
f1  <- H.mix.bs(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7))
f2  <- H.mix.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7),nu = 4, mfamily = "t")
f3  <- H.mix.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7),nu = 4, mfamily = "sl")
f4  <- H.mix.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
#legend("topleft",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("hazardgdif3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,13,0.01)
f1  <- H.mix.bs(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6))
f2  <- H.mix.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6),nu = 4, mfamily = "t")
f3  <- H.mix.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6),nu = 4, mfamily = "sl")
f4  <- H.mix.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
#legend("topleft",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("hazardgdif4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,8,0.01)
f1  <- H.mix.bs(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3))
f2  <- H.mix.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3),nu = 4, mfamily = "t")
f3  <- H.mix.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3),nu = 4, mfamily = "sl")
f4  <- H.mix.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
#legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

