####################################################################################################################
#----------------------------------        Figure FM-BS-SMN      --------------------------------------------------#
####################################################################################################################
#Estamos considerandos tres escenarios
#Escenario 1: Mal separados
#Escenario 2: Medianamente separados
#Escenario 3: Bien separados
####################################################################################################################
#Cargando pacotes
setwd("C:/Users/Luiz/Dropbox/Pacotes/bssn/R")

source("functions.bssmn.R")
source("functions.mixbssmn.R")
source("functions.mixbs.R")
source("functions.R")
####################################################################################################################

postscript("dens_escenario1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,3,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.6,0.4), alpha=c(0.75,0.5), beta=c(0.5,1.5))
f2  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.75,0.5), beta=c(0.5,1.5),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.75,0.5), beta=c(0.5,1.5),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.75,0.5), beta=c(0.5,1.5),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript


postscript("dens_escenario2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,15,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.6,0.4), alpha=c(0.5,0.5), beta=c(1,6))
f2  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.5), beta=c(1,6),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.5), beta=c(1,6),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.5), beta=c(1,6),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript


postscript("dens_escenario3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,8,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,5))
f2  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,5),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,5),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,5),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript




#########################################################################
postscript("densgdif1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,30,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7))
f2  <- d.mixed.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.5,0.75), beta=c(3,7),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("densgdif2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,30,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7))
f2  <- d.mixed.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.2,0.8), alpha=c(0.25,0.35), beta=c(3,7),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("densgdif3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,10,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6))
f2  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.5,0.25), beta=c(1,6),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("densgdif4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,10,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3))
f2  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.25,0.25), beta=c(1,3),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("densgdif5.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,10,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.6,0.4), alpha=c(1.25,0.25), beta=c(3,1))
f2  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(1.25,0.25), beta=c(3,1),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(1.25,0.25), beta=c(3,1),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(1.25,0.25), beta=c(3,1),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("densgdif6.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(0.01,30,0.001)
f1  <- d.mixed.bs(ti,pii=c(0.6,0.4), alpha=c(0.75,0.25), beta=c(10,3))
f2  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.75,0.25), beta=c(10,3),nu = 4, mfamily = "t")
f3  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.75,0.25), beta=c(10,3),nu = 4, mfamily = "sl")
f4  <- d.mixed.bssmn(ti,pii=c(0.6,0.4), alpha=c(0.75,0.25), beta=c(10,3),nu = c(0.1,0.1), mfamily = "cn")
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend("topright",c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript


