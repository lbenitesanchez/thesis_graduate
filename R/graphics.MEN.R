#install.packages("VGAM")
#library(VGAM)#tem a funcao dslash

dslash <- function(x, mu = 0, sigma = 1, log = FALSE,
                   smallno = .Machine$double.eps * 1000) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  #if (!is.Numeric(sigma) || any(sigma <= 0))
   # stop("argument 'sigma' must be positive")
  L <- max(length(x), length(mu), length(sigma))
  if (length(x)     != L) x     <- rep_len(x,     L)
  if (length(mu)    != L) mu    <- rep_len(mu,    L)
  if (length(sigma) != L) sigma <- rep_len(sigma, L)

  zedd <- (x-mu)/sigma
  if (log.arg) {
    ifelse(abs(zedd) < smallno,
           -log(2*sigma*sqrt(2*pi)),
           log1p(-exp(-zedd^2/2)) - log(sqrt(2*pi)*sigma*zedd^2))
  } else {
    ifelse(abs(zedd) < smallno,
           1/(2*sigma*sqrt(2*pi)),
           -expm1(-zedd^2/2)/(sqrt(2*pi)*sigma*zedd^2))
  }
}

dcn <- function(x, mu = 0, sigma =1, nu = 0.1, gamma=0.1)
{
 return(nu*dnorm(x,mu,gamma^{-1}*sigma) + (1-nu)*dnorm(x,mu,sigma))

}

postscript("densInifig2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(-35,20,0.001)
f1  <- dnorm(ti,1,2)
f2  <- dt(ti,df = 3)
f3  <- dslash(ti,1,2)
f4  <- dcn(ti,1,2,nu=0.1,gamma=0.1)
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend(-36,0.35,c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript

postscript("densInifig1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ti  <- seq(-27,10,0.001)
f1  <- dnorm(ti)
f2  <- dt(ti,df = 4)
f3  <- dslash(ti)
f4  <- dcn(ti)
den <- cbind(f1,f2,f3,f4)
matplot(ti,den,type = "l", col = c("black","red","blue","green"), ylab = "f(t)" , xlab="t",lwd=2,lty=c(1,2,3,4))
legend(-28,0.35,c("Normal","T","Slash","Normal contaminada"),col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=2,bty="n",inset=0.1,cex=1.0, pt.cex = 1)
dev.off() #Fechando o dispositivo potscript
