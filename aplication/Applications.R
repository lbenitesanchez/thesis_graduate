library(sn)
setwd("~/Dropbox/Pacotes/bssn/R")
setwd("C:/Users/Luiz/Dropbox/Pacotes/bssn/R")
setwd("C:/Users/Larissa/Dropbox/Pacotes/bssn/R")
setwd("C:/Users/luis_/Dropbox/Luis y Rocio/Research/Pacotes/bssn/R")


source("functions.R")
source("functions.mixbssmn.R")
source("functions.bssmn.R")
source("functions.mixbs.R")

source("algEMmixbssmn.R")
source("algEMmixbs.R")
source("algEMmixbs.R")
library(bssn)

#Enzyme  #Dados utilizados por Leiva
###############################################
xi = ti = c(0.021, 0.031, 0.044, 0.061, 0.070, 0.077, 0.078, 0.080, 0.081, rep(0.083,2), 0.084, 0.085, 0.088, rep(0.096,2), 0.100, 0.106, rep(0.108,2), 0.109, 0.111, 0.112, 0.113, rep(0.118,2), 0.122,
            rep(0.124,3), rep(0.126,3), 0.128, 0.129, rep(0.130,2), 0.131, rep(0.132,2), 0.134, rep(0.137,2), 0.138, rep(0.142,2), 0.144, rep(0.148,2), rep(0.149,2), 0.151,
            rep(0.152,2), 0.157, 0.159, rep(0.162,3), rep(0.166,2), rep(0.167,2), 0.171, rep(0.172,2), 0.174, 0.175, rep(0.176,2), 0.178, rep(0.179,2), rep(0.180,2), 0.182, 0.183,
            rep(0.184,2), rep(0.185,2), rep(0.190,2), 0.191, rep(0.192,4), 0.193, 0.194, 0.195, rep(0.198,3), rep(0.200,2),0.204, 0.205, 0.207, 0.210, 0.213, 0.214, 0.215,
            0.216, 0.222, 0.224, 0.225, rep(0.230,2), rep(0.232,2), 0.236, 0.238, 0.240, 0.241, 0.246, 0.250, 0.255, rep(0.258,3), 0.263, 0.264, 0.265, 0.275,
            0.277, 0.279, 0.280, rep(0.291,2), 0.292, 0.293, 0.299, rep(0.305,2), 0.308, 0.309, 0.313, rep(0.320,2), 0.325, 0.340, 0.342, 0.347, 0.357, 0.360,
            0.368, 0.379, 0.387, 0.408, 0.409, 0.466, 0.520, 0.673, 0.687, 0.709, 0.718, 0.818, 0.839, 0.853, 0.860, 0.867, 0.875, 0.895, 0.902,
            0.909, 0.923, 0.929, 0.933, 0.945, rep(0.953,2), 0.962, 0.967, 0.978, 0.985, 0.998, 1.003, 1.004, 1.007, 1.010, 1.018, rep(1.036,2), 1.040, 1.052,
            1.073, 1.075, 1.113, 1.123, 1.161, 1.166, 1.169, 1.173, 1.176, rep(1.195,2), 1.200, 1.229, 1.231, 1.241, 1.261, 1.271, 1.279, 1.293, 1.298,
            1.300, 1.310, 1.326, 1.329, 1.361, 1.369, 1.371, 1.374, 1.405, 1.419, 1.452, 1.461, 1.488, 1.516, 1.567, 1.573, 1.604, 1.622, 1.625,
            1.633, 1.672, 1.713, 1.741, 1.742, 1.768, 1.811, 1.848, 1.861, 1.885, 1.944, 1.950, 2.016, 2.183, 2.264, 2.338, 2.427, 2.518, 2.545,2.880)

postscript("densityENZYME.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(density(ti),main="", xlab="")
dev.off() #Fechando o dispositivo potscript

postscript("boxplotENZYME.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
boxplot(ti)
dev.off() #Fechando o dispositivo potscript

fitBHg1     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=1, get.init = "BH", accuracy = 10^-6, iter.max = 100, kmeans.param = NULL,aitken=FALSE, rates=FALSE)
fitBHg2     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=2, get.init = "BH", accuracy = 10^-6, iter.max = 100, kmeans.param = NULL,aitken=FALSE, rates=FALSE)
fitBHg3     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=3, get.init = "BH", accuracy = 10^-6, iter.max = 100, kmeans.param = NULL,aitken=FALSE, rates=FALSE)
fitBHg4     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=4, get.init = "BH", accuracy = 10^-6, iter.max = 100, kmeans.param = NULL,aitken=FALSE, rates=FALSE)

initialTg1  <- initial.Values(ti, g=1, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)
initialTg2  <- initial.Values(ti, g=2, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = TRUE)
initialTg3  <- initial.Values(ti, g=3, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)
initialTg4  <- initial.Values(ti, g=4, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)

fit_t_g1    <- EMmixbssmn2(ti,initialTg1$alpha, initialTg1$beta, initialTg1$nu, initialTg1$pii, g=1, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g2    <- EMmixbssmn2(ti,initialTg2$alpha, initialTg2$beta, initialTg2$nu, initialTg2$pii, g=2, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g3    <- EMmixbssmn2(ti,initialTg3$alpha, initialTg3$beta, initialTg3$nu, initialTg3$pii, g=3, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g4    <- EMmixbssmn2(ti,initialTg4$alpha, initialTg4$beta, initialTg4$nu, initialTg4$pii, g=4, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)

round(c(fit_t_g1$result$lk   ,fit_t_g2$result$lk   ,fit_t_g3$result$lk   , fit_t_g4$result$lk),digits = 4)
c(fit_t_g1$result$aic  ,fit_t_g2$result$aic  ,fit_t_g3$result$aic  , fit_t_g4$result$aic)
c(fit_t_g1$result$bic  ,fit_t_g2$result$bic  ,fit_t_g3$result$bic  , fit_t_g4$result$bic)
c(fit_t_g1$result$iter,fit_t_g2$result$iter,fit_t_g3$result$iter, fit_t_g4$result$iter)
round(c(fit_t_g1$result$rates,fit_t_g2$result$rates,fit_t_g3$result$rates, fit_t_g4$result$rates),digits = 4)

round(InfmatrixmixBSSMN1(ti,fit_t_g2$result$pii,fit_t_g2$result$alpha,fit_t_g2$result$beta,fit_t_g2$result$nu,mfamily="t")$EP,digits = 4)

#-----------------------------------------------------------------------------------------------------
initialSLg1  <- initial.Values(ti,g=1,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialSLg2  <- initial.Values(ti,g=2,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = TRUE)
initialSLg3  <- initial.Values(ti,g=3,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialSLg4  <- initial.Values(ti,g=4,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)

fit_sl_g1  <- EMmixbssmn2(ti,initialSLg1$alpha, initialSLg1$beta, initialSLg1$nu, initialSLg1$pii, g=1, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE,rates=FALSE)
fit_sl_g2  <- EMmixbssmn2(ti,initialSLg2$alpha, initialSLg2$beta, initialSLg2$nu, initialSLg2$pii, g=2, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_sl_g3  <- EMmixbssmn2(ti,initialSLg3$alpha, initialSLg3$beta, initialSLg3$nu, initialSLg3$pii, g=3, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_sl_g4  <- EMmixbssmn2(ti,initialSLg4$alpha, initialSLg4$beta, initialSLg4$nu, initialSLg4$pii, g=4, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE)

c(fit_sl_g1$result$lk,fit_sl_g2$result$lk,fit_sl_g3$result$lk,fit_sl_g4$result$lk)
c(fit_sl_g1$result$aic,fit_sl_g2$result$aic,fit_sl_g3$result$aic,fit_sl_g4$result$aic)
c(fit_sl_g1$result$bic,fit_sl_g2$result$bic,fit_sl_g3$result$bic,fit_sl_g4$result$bic)
c(fit_sl_g1$result$iter,fit_sl_g2$result$iter,fit_sl_g3$result$iter,fit_sl_g4$result$iter)
c(fit_sl_g1$result$rates,fit_sl_g2$result$rates,fit_sl_g3$result$rates,fit_sl_g4$result$rates)

round(InfmatrixmixBSSMN1(ti,fit_sl_g2$result$pii,fit_sl_g2$result$alpha,fit_sl_g2$result$beta,fit_sl_g2$result$nu,mfamily="t")$EP,digits = 4)


#-----------------------------------------------------------------------------------------------------
initialCNg1  <- initial.Values(ti,g=1,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialCNg2  <- initial.Values(ti,g=2,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = TRUE)
initialCNg3  <- initial.Values(ti,g=3,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialCNg4  <- initial.Values(ti,g=4,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)

fit_cn_g1  <- EMmixbssmn2(ti,initialCNg1$alpha, initialCNg1$beta, initialCNg1$nu, initialCNg1$pii, g=1, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE,rates=FALSE)
fit_cn_g2  <- EMmixbssmn2(ti,initialCNg2$alpha, initialCNg2$beta, initialCNg2$nu, initialCNg2$pii, g=2, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_cn_g3  <- EMmixbssmn2(ti,initialCNg3$alpha, initialCNg3$beta, initialCNg3$nu, initialCNg3$pii, g=3, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_cn_g4  <- EMmixbssmn2(ti,initialCNg4$alpha, initialCNg4$beta, initialCNg4$nu, initialCNg4$pii, g=4, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE)

InfmatrixmixBSSMN1(ti,fit_cn_g2$result$pii,fit_cn_g2$result$alpha,fit_cn_g2$result$beta,fit_cn_g2$result$nu,mfamily="cn")

c(fit_cn_g1$result$lk,fit_cn_g2$result$lk,fit_cn_g3$result$lk,fit_cn_g4$result$lk)
c(fit_cn_g1$result$aic,fit_cn_g2$result$aic,fit_cn_g3$result$aic,fit_cn_g4$result$aic)
c(fit_cn_g1$result$bic,fit_cn_g2$result$bic,fit_cn_g3$result$bic,fit_cn_g4$result$bic)
c(fit_cn_g1$result$iter,fit_cn_g2$result$iter,fit_cn_g3$result$iter,fit_cn_g4$result$iter)
c(fit_cn_g1$result$rates,fit_cn_g2$result$rates,fit_cn_g3$result$rates,fit_sl_g4$result$rates)

######################################################################################################
#                                         Enzyme dataset
######################################################################################################
postscript("MIX.FM.BS.SMN.densityENZYMEg=1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, main="g=1",xlab="t",col = "grey", border=FALSE,ylim=c(0,4))
xx  <- seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1))
lines(xx, d.mixed.bs(xx, fitBHg1$result$pii, fitBHg1$result$alpha, fitBHg1$result$beta), col="red", lty=2,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g1$result$pii, fit_t_g1$result$alpha, fit_t_g1$result$beta, fit_t_g1$result$nu,mfamily="t"), col="blue", lty=2,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g1$result$pii, fit_cn_g1$result$alpha, fit_cn_g1$result$beta, fit_cn_g1$result$nu,mfamily="cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g1$result$pii, fit_sl_g1$result$alpha, fit_sl_g1$result$beta, fit_sl_g1$result$nu,mfamily="sl"), col="green", lty=2,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript


postscript("MIX.FM.BS.SMN.densityENZYMEg=2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, main="g=2",xlab="t",col = "grey", border=FALSE,ylim=c(0,4))
xx  <- seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1))
lines(xx, d.mixed.bs(xx, fitBHg2$result$pii, fitBHg2$result$alpha, fitBHg2$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g2$result$pii, fit_t_g2$result$alpha, fit_t_g2$result$beta, fit_t_g2$result$nu,mfamily="t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g2$result$pii, fit_cn_g2$result$alpha, fit_cn_g2$result$beta, fit_cn_g2$result$nu,mfamily="cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g2$result$pii, fit_sl_g2$result$alpha, fit_sl_g2$result$beta, fit_sl_g2$result$nu,mfamily="sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black","green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityENZYMEg=3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, main="g=3",xlab="t",col = "grey", border=FALSE,ylim=c(0,4))
xx  <- seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1))
lines(xx, d.mixed.bs(xx, fitBHg3$result$pii, fitBHg3$result$alpha, fitBHg3$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g3$result$pii, fit_t_g3$result$alpha, fit_t_g3$result$beta, fit_t_g3$result$nu,mfamily="t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g3$result$pii, fit_sl_g3$result$alpha, fit_sl_g3$result$beta, fit_sl_g3$result$nu,mfamily="sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","SL"),col=c("red", "blue","green"),lty=c(2,2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityENZYMEg=4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, main="g=4",xlab="t",col = "grey", border=FALSE,ylim=c(0,4))
xx  <- seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1))
lines(xx, d.mixed.bs(xx, fitBHg4$result$pii, fitBHg4$result$alpha, fitBHg4$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g4$result$pii, fit_t_g4$result$alpha, fit_t_g4$result$beta, fit_t_g4$result$nu,mfamily="t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g4$result$pii, fit_sl_g4$result$alpha, fit_sl_g4$result$beta, fit_sl_g4$result$nu,mfamily="sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","SL"),col=c("red", "blue","green"),lty=c(2,2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

################################################33
#Acumulada
postscript("cdfsimmixBMI_g1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=1",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,10,0.01)
df1    <- p.mixed.bs(xx,pii=fitBHg1$result$pii,alpha=fitBHg1$result$alpha,beta=fitBHg1$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g1$result$pii,alpha=fit_t_g1$result$alpha,beta=fit_t_g1$result$beta,nu = fit_t_g1$result$nu,mfamily="t")
df3    <- p.mixed.bssmn(xx,pii=fit_cn_g1$result$pii,alpha=fit_cn_g1$result$alpha,beta=fit_cn_g1$result$beta,nu = fit_cn_g1$result$nu,mfamily="cn")
df4    <- p.mixed.bssmn(xx,pii=fit_sl_g1$result$pii,alpha=fit_sl_g1$result$alpha,beta=fit_sl_g1$result$beta,nu = fit_sl_g1$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
lines(xx, df3,col="black", lty=1,lwd=2)
lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("cdfsimmixBMI_g2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=2",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,10,0.01)
df1    <- p.mixed.bs(xx,pii=fitBHg2$result$pii,alpha=fitBHg2$result$alpha,beta=fitBHg2$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g2$result$pii,alpha=fit_t_g2$result$alpha,beta=fit_t_g2$result$beta,nu = fit_t_g2$result$nu,mfamily="t")
df3    <- p.mixed.bssmn(xx,pii=fit_cn_g2$result$pii,alpha=fit_cn_g2$result$alpha,beta=fit_cn_g2$result$beta,nu = fit_cn_g2$result$nu,mfamily="cn")
#df4    <- p.mixed.bssmn(xx,pii=fit_sl_g2$result$pii,alpha=fit_sl_g2$result$alpha,beta=fit_sl_g2$result$beta,nu = fit_sl_g2$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
lines(xx, df3,col="black", lty=1,lwd=2)
#lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T","NC"),col=c("red", "blue","black"),lty=c(2,2,1),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("cdfsimmixBMI_g3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=3",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,20,0.01)
df1    <- p.mixed.bs(xx,pii=fitBHg3$result$pii,alpha=fitBHg3$result$alpha,beta=fitBHg3$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g3$result$pii,alpha=fit_t_g3$result$alpha,beta=fit_t_g3$result$beta,nu = fit_t_g2$result$nu,mfamily="t")
#df4    <- p.mixed.bssmn(xx,pii=fit_sl_g3$result$pii,alpha=fit_sl_g3$result$alpha,beta=fit_sl_g3$result$beta,nu = fit_sl_g2$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
#lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T"),col=c("red", "blue"),lty=c(2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("cdfsimmixBMI_g4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=4",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,20,0.01)
df1    <- p.mixed.bs(xx,pii=fitBHg4$result$pii,alpha=fitBHg4$result$alpha,beta=fitBHg4$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g4$result$pii,alpha=fit_t_g4$result$alpha,beta=fit_t_g4$result$beta,nu = fit_t_g2$result$nu,mfamily="t")
#df4    <- p.mixed.bssmn(xx,pii=fit_sl_g4$result$pii,alpha=fit_sl_g4$result$alpha,beta=fit_sl_g4$result$beta,nu = fit_sl_g2$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
#lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T"),col=c("red", "blue"),lty=c(2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript


##########################################################################
#Survival
#survival
library(survival)

fit     <- survfit(Surv(ti) ~ 1)
postscript("survival_apliBMIg1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=1")
xx     <- seq(0.1,length(ti),0.1)
df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg1$result$pii,alpha=fitBHg1$result$alpha,beta=fitBHg1$result$beta)
df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g1$result$pii,alpha=fit_t_g1$result$alpha,beta=fit_t_g1$result$beta,nu=fit_t_g1$result$nu,mfamily = "t")
df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g1$result$pii,alpha=fit_cn_g1$result$alpha,beta=fit_cn_g1$result$beta,nu=fit_cn_g1$result$nu,mfamily = "cn")
df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g1$result$pii,alpha=fit_sl_g1$result$alpha,beta=fit_sl_g1$result$beta,nu=fit_sl_g1$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
lines(xx, df_g_3,col="black", lty=1,lwd=2)
lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript


postscript("survival_apliBMIg2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=2")
xx     <- seq(0.1,length(ti),0.1)
df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg2$result$pii,alpha=fitBHg2$result$alpha,beta=fitBHg2$result$beta)
df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g2$result$pii,alpha=fit_t_g2$result$alpha,beta=fit_t_g2$result$beta,nu=fit_t_g2$result$nu,mfamily = "t")
#df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g2$result$pii,alpha=fit_cn_g2$result$alpha,beta=fit_cn_g2$result$beta,nu=fit_cn_g2$result$nu,mfamily = "cn")
#df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g2$result$pii,alpha=fit_sl_g2$result$alpha,beta=fit_sl_g2$result$beta,nu=fit_sl_g2$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
#lines(xx, df_g_3,col="black", lty=1,lwd=2)
#lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T"),col=c("red", "blue"),lty=c(2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("survival_apliBMIg3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=3")
xx     <- seq(0.1,length(ti),0.1)
df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg3$result$pii,alpha=fitBHg3$result$alpha,beta=fitBHg3$result$beta)
df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g3$result$pii,alpha=fit_t_g3$result$alpha,beta=fit_t_g3$result$beta,nu=fit_t_g3$result$nu,mfamily = "t")
#df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g3$result$pii,alpha=fit_cn_g3$result$alpha,beta=fit_cn_g3$result$beta,nu=fit_cn_g3$result$nu,mfamily = "cn")
#df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g3$result$pii,alpha=fit_sl_g3$result$alpha,beta=fit_sl_g3$result$beta,nu=fit_sl_g3$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
#lines(xx, df_g_3,col="black", lty=1,lwd=2)
#lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T"),col=c("red", "blue"),lty=c(2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("survival_apliBMIg4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=4")
xx     <- seq(0.1,length(ti),0.1)
df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg4$result$pii,alpha=fitBHg4$result$alpha,beta=fitBHg4$result$beta)
df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g4$result$pii,alpha=fit_t_g4$result$alpha,beta=fit_t_g4$result$beta,nu=fit_t_g4$result$nu,mfamily = "t")
#df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g3$result$pii,alpha=fit_cn_g3$result$alpha,beta=fit_cn_g3$result$beta,nu=fit_cn_g3$result$nu,mfamily = "cn")
#df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g3$result$pii,alpha=fit_sl_g3$result$alpha,beta=fit_sl_g3$result$beta,nu=fit_sl_g3$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
#lines(xx, df_g_3,col="black", lty=1,lwd=2)
#lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T"),col=c("red", "blue"),lty=c(2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

##############################################################################
#Tabla de Estimativas
Estimates <- round(c(fit_cn_g2$result$alpha, fit_cn_g2$result$beta, fit_cn_g2$result$pii[1]),digits=4)
EP        <- round(InfmatrixmixBSSMN1(ti,fit_cn_g2$result$pii,fit_cn_g2$result$alpha,fit_cn_g2$result$beta,fit_cn_g2$result$nu,mfamily="cn")$EP,digits=4)
U         <- round(Estimates + 1.96*EP,digits=4)
L         <- round(Estimates - 1.96*EP,digits=4)


Estimates <- round(c(fit_t_g2$result$alpha, fit_t_g2$result$beta, fit_t_g2$result$pii[2]),digits=4)
EP        <- round(InfmatrixmixBSSMN1(ti,fit_t_g2$result$pii,fit_t_g2$result$alpha,fit_t_g2$result$beta,fit_t_g2$result$nu,mfamily="t")$EP,digits=4)
U         <- round(Estimates + 1.96*EP,digits=4)
L         <- round(Estimates - 1.96*EP,digits=4)

Estimates <- round(c(fit_sl_g2$result$alpha, fit_sl_g2$result$beta, fit_sl_g2$result$pii[2]),digits=4)
EP        <- round(InfmatrixmixBSSMN1(ti,fit_sl_g2$result$pii,fit_sl_g2$result$alpha,fit_sl_g2$result$beta,fit_sl_g2$result$nu,mfamily="sl")$EP,digits=4)
U         <- round(Estimates + 1.96*EP,digits=4)
L         <- round(Estimates - 1.96*EP,digits=4)


#################################################################################3
#
#             Figuras de las densidades con sus respectivas flechas
#
#################################################################################
postscript("MIX.FM.BS.SMN.densityENZYMEg=1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,3.50),main="")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg1$result$pii, fitBHg1$result$alpha, fitBHg1$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g1$result$pii, fit_t_g1$result$alpha, fit_t_g1$result$beta, fit_t_g1$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g1$result$pii, fit_cn_g1$result$alpha, fit_cn_g1$result$beta, fit_cn_g1$result$nu, mfamily = "cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g1$result$pii, fit_sl_g1$result$alpha, fit_sl_g1$result$beta, fit_sl_g1$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
arrows(0.6201, 2.1362, 0.2005, 1.579, xpd = TRUE,length=0.1, col="red", lty=1, lwd=2)
arrows(0.3647, 2.7682, 0.1134, 1.901, xpd = TRUE,length=0.1, col="blue", lty=1, lwd=2)
arrows(0.7416, 1.3680, 0.3321, 1.094, xpd = TRUE,length=0.1, col="black", lty=1, lwd=2)
arrows(0.9511, 0.6983, 0.6625, 0.500, xpd = TRUE,length=0.1, col="green", lty=1, lwd=2)
text(1.45, 2.1990, expression(paste(alpha[1],"=1.1458",", ",beta[1],"=0.3783")), pos = 2,cex=0.8, col="red")
text(1.49, 2.8291, expression(paste(alpha[1],"=1.1115",", ",beta[1],"=0.3788",", ",nu,"=20")), pos = 2,cex=0.8, col="blue")
text(2.25, 1.3984, expression(paste(alpha[1],"=1.1396",", ",beta[1],"=0.3783",", ",nu,"=0.9900",", ",gamma,"=0.9892")), pos = 2,cex=0.8, col="black")
text(1.97, 0.7287, expression(paste(alpha[1],"=1.1171",", ",beta[1],"=0.3784",", ",nu,"=20")), pos = 2,cex=0.8, col="green")
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityENZYMEg=2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,3.50),main="")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg2$result$pii, fitBHg2$result$alpha, fitBHg2$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g2$result$pii, fit_t_g2$result$alpha, fit_t_g2$result$beta, fit_t_g2$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g2$result$pii, fit_cn_g2$result$alpha, fit_cn_g2$result$beta, fit_cn_g2$result$nu, mfamily = "cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g2$result$pii, fit_sl_g2$result$alpha, fit_sl_g2$result$beta, fit_sl_g2$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
arrows(0.6102, 2.1138, 0.2165, 1.977, xpd = TRUE,length=0.1, col="red", lty=1, lwd=2)
arrows(0.3647, 2.7682, 0.1134, 1.901, xpd = TRUE,length=0.1, col="blue", lty=1, lwd=2)
arrows(0.7414, 1.3528, 0.2907, 1.064, xpd = TRUE,length=0.1, col="black", lty=1, lwd=2)
arrows(0.6729, 0.6374, 0.4675, 0.150, xpd = TRUE,length=0.1, col="green", lty=1, lwd=2)
text(2.09, 2.1590, expression(paste(p[1],"=0.6259",", ",alpha[1],"=0.5238",", ",beta[1],"=0.1734")), pos = 2,cex=0.8, col="red")
text(2.09, 1.9590, expression(paste(alpha[2],"=0.3231",", ",beta[2],"=1.2669")), pos = 2,cex=0.8, col="red")
text(1.75, 2.8891, expression(paste(p[1],"=0.6267",alpha[1],"=0.4085",", ",beta[1],"=0.3788")), pos = 2,cex=0.8, col="blue")
text(1.75, 2.6591, expression(paste(alpha[2],"=0.2859",", ",beta[2],"=1.2669",", ",nu,"=20")), pos = 2,cex=0.8, col="blue")
text(2.66, 1.3984, expression(paste(p[1],"=0.6263",", ",alpha[1],"=0.4263",", ",beta[1],"=0.1812")), pos = 2,cex=0.8, col="black")
text(2.66, 1.1784, expression(paste(alpha[2],"=0.3123",", ",beta[2],"=1.2573",", ",nu,"=0.0555",", ",gamma,"=0.0861")), pos = 2,cex=0.8, col="black")
text(2.52, 0.8287, expression(paste(p[1],"=0.6257",", ",alpha[1],"=0.3561",", ",beta[1],"=0.1806",", ",nu,"=2.01")), pos = 2,cex=0.8, col="green")
text(2.52, 0.6187, expression(paste(alpha[2],"=0.2512",", ",beta[2],"=1.2523")), pos = 2,cex=0.8, col="green")
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityENZYMEg=3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,3.50),main="")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg3$result$pii, fitBHg3$result$alpha, fitBHg3$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g3$result$pii, fit_t_g3$result$alpha, fit_t_g3$result$beta, fit_t_g3$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g3$result$pii, fit_sl_g3$result$alpha, fit_sl_g3$result$beta, fit_sl_g3$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","SL"),col=c("red", "blue", "green"),lty=c(2,2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
arrows(0.6102, 2.1138, 0.2165, 1.977, xpd = TRUE,length=0.1, col="red", lty=1, lwd=2)
arrows(0.3647, 2.7682, 0.1134, 1.901, xpd = TRUE,length=0.1, col="blue", lty=1, lwd=2)
arrows(0.6729, 0.6374, 0.4675, 0.150, xpd = TRUE,length=0.1, col="green", lty=1, lwd=2)
text(2.09, 2.1590, expression(paste(p[1],"=0.6277",", ",alpha[1],"=0.5281",", ",beta[1],"=0.1742")), pos = 2,cex=0.8, col="red")
text(2.09, 1.9590, expression(paste(p[2],"=0.2605",", ",alpha[2],"=0.2130",", ",beta[2],"=1.0921")), pos = 2,cex=0.8, col="red")
text(2.09, 1.7590, expression(paste(alpha[3],"=0.2283",", ",beta[3],"=1.8020")), pos = 2,cex=0.8, col="red")

text(1.95, 3.1, expression(paste(p[1],"=0.1328",", ",alpha[1],"=0.2230",", ",beta[1],"=1.6735")), pos = 2,cex=0.8, col="blue")
text(1.95, 2.88, expression(paste(p[2],"=0.2368",", ",alpha[2],"=0.1848",", ",beta[2],"=1.0827")), pos = 2,cex=0.8, col="blue")
text(1.95, 2.66, expression(paste(alpha[3],"=0.4095",", ",beta[3],"=0.1766",", ",nu,"=4.9509")), pos = 2,cex=0.8, col="blue")

text(2.52, 1.0287, expression(paste(p[1],"=0.1224",", ",alpha[1],"=0.1882",", ",beta[1],"=1.7308",", ",nu,"=2.0100")), pos = 2,cex=0.8, col="green")
text(2.52, 0.8187, expression(paste(p[2],"=0.6282",", ",alpha[2],"=0.3599",", ",beta[2],"=0.1753")), pos = 2,cex=0.8, col="green")
text(2.52, 0.6187, expression(paste(alpha[3],"=0.1642",", ",beta[3],"=1.0840")), pos = 2,cex=0.8, col="green")
dev.off() #Fechando o dispositivo potscript


postscript("MIX.FM.BS.SMN.densityENZYMEg=4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,3.50),main="")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg4$result$pii, fitBHg4$result$alpha, fitBHg4$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g4$result$pii, fit_t_g4$result$alpha, fit_t_g4$result$beta, fit_t_g4$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g4$result$pii, fit_sl_g4$result$alpha, fit_sl_g4$result$beta, fit_sl_g4$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","SL"),col=c("red", "blue", "green"),lty=c(2,2,2),lwd=2,bty="n",inset=0.1,cex=0.8)
arrows(0.6102, 2.1138, 0.2165, 1.977, xpd = TRUE,length=0.1, col="red", lty=1, lwd=2)
arrows(0.3647, 2.7682, 0.1134, 1.901, xpd = TRUE,length=0.1, col="blue", lty=1, lwd=2)
arrows(0.6729, 0.6374, 0.4675, 0.150, xpd = TRUE,length=0.1, col="green", lty=1, lwd=2)
text(2.09, 2.1590, expression(paste(p[1],"=0.6277",", ",alpha[1],"=0.5281",", ",beta[1],"=0.1741")), pos = 2,cex=0.8, col="red")
text(2.09, 1.9590, expression(paste(p[2],"=0.2925",", ",alpha[2],"=0.2236",", ",beta[2],"=1.1238")), pos = 2,cex=0.8, col="red")
text(2.09, 1.7590, expression(paste(p[3],"=0.0513",", ",alpha[3],"=0.0819",", ",beta[3],"=1.7713")), pos = 2,cex=0.8, col="red")
text(2.09, 1.5590, expression(paste(alpha[4],"=0.0923",", ",beta[4],"=2.4318")), pos = 2,cex=0.8, col="red")

text(1.95, 3.32, expression(paste(p[1],"=0.1328",", ",alpha[1],"=0.2230",", ",beta[1],"=1.6735")), pos = 2,cex=0.8, col="blue")
text(1.95, 3.10, expression(paste(p[2],"=0.2368",", ",alpha[2],"=0.1848",", ",beta[2],"=1.0827")), pos = 2,cex=0.8, col="blue")
text(1.95, 2.88, expression(paste(p[3],"=0.2368",", ",alpha[3],"=0.1848",", ",beta[3],"=1.0827")), pos = 2,cex=0.8, col="blue")
text(1.95, 2.66, expression(paste(alpha[4],"=0.4095",", ",beta[4],"=0.1766",", ",nu,"=4.9509")), pos = 2,cex=0.8, col="blue")

text(2.52, 1.2287, expression(paste(p[1],"=0.5044",", ",alpha[1],"=0.3428",", ",beta[1],"=0.1589",", ",nu,"=2.0100")), pos = 2,cex=0.8, col="green")
text(2.52, 1.0187, expression(paste(p[2],"=0.2715",", ",alpha[2],"=0.1744",", ",beta[2],"=1.1017")), pos = 2,cex=0.8, col="green")
text(2.52, 0.8187, expression(paste(p[3],"=0.1027",", ",alpha[3],"=0.1768",", ",beta[3],"=1.8025")), pos = 2,cex=0.8, col="green")
text(2.52, 0.6187, expression(paste(alpha[4],"=0.2245",", ",beta[4],"=0.2599")), pos = 2,cex=0.8, col="green")
dev.off() #Fechando o dispositivo potscript
#################################################################################


######################################################################################################
#                                         BMI dataset
######################################################################################################
library(mixsmsn)
data("bmi")
ti<- bmi[,1]
hist(ti)
plot(density(ti))

postscript("densityBMI.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(density(ti),main="", xlab="")
dev.off() #Fechando o dispositivo potscript

postscript("boxplotBMI.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
boxplot(ti)
dev.off() #Fechando o dispositivo potscript


fitBHg1     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=1, get.init = "BH", accuracy = 10^-6, iter.max = 100, kmeans.param = NULL,aitken=FALSE, rates=TRUE)
fitBHg2     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=2, get.init = "BH", accuracy = 10^-6, iter.max = 100, kmeans.param = NULL,aitken=FALSE, rates=TRUE)
fitBHg3     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=3, get.init = "BH", accuracy = 10^-6, iter.max = 500, kmeans.param = NULL,aitken=FALSE, rates=TRUE)
fitBHg4     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=4, get.init = "BH", accuracy = 10^-6, iter.max = 500, kmeans.param = NULL,aitken=FALSE, rates=TRUE)
fitBHg5     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=5, get.init = "BH", accuracy = 10^-6, iter.max = 500, kmeans.param = NULL,aitken=FALSE, rates=TRUE)
fitBHg6     <- alg.EM.mix.bs(ti,  alpha = NULL, beta = NULL, pii = NULL, g=6, get.init = "BH", accuracy = 10^-6, iter.max = 500, kmeans.param = NULL,aitken=FALSE, rates=TRUE)


round(c(fitBHg1$result$lk   ,fitBHg2$result$lk   ,fitBHg3$result$lk   , fitBHg4$result$lk),digits = 4)
c(fitBHg1$result$aic  ,fitBHg2$result$aic  ,fitBHg3$result$aic  , fitBHg4$result$aic)
c(fitBHg1$result$bic  ,fitBHg2$result$bic  ,fitBHg3$result$bic  , fitBHg4$result$bic)
c(fitBHg1$result$iter,fitBHg2$result$iter,fitBHg3$result$iter, fitBHg4$result$iter)
round(c(fitBHg1$result$rates,fitBHg2$result$rates,fitBHg3$result$rates, fitBHg4$result$rates),digits = 4)


round(InfmatrixmixBS2(ti,fitBHg2$result$pii,fitBHg2$result$alpha,fitBHg2$result$beta)$EP,digits = 4)


initialTg1  <- initial.Values(ti, g=1, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)
initialTg2  <- initial.Values(ti, g=2, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = TRUE)
initialTg3  <- initial.Values(ti, g=3, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)
initialTg4  <- initial.Values(ti, g=4, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)
initialTg5  <- initial.Values(ti, g=5, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)
initialTg6  <- initial.Values(ti, g=6, algorithm="k-means", mfamily="t", lower=1, upper=20, space=0.01, plotLog = TRUE, searchNU=TRUE, printNU=FALSE, saveFigure = FALSE)


fit_t_g1    <- EMmixbssmn2(ti,initialTg1$alpha, initialTg1$beta, initialTg1$nu, initialTg1$pii, g=1, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g2    <- EMmixbssmn2(ti,initialTg2$alpha, initialTg2$beta, initialTg2$nu, initialTg2$pii, g=2, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g3    <- EMmixbssmn2(ti,initialTg3$alpha, initialTg3$beta, initialTg3$nu, initialTg3$pii, g=3, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g4    <- EMmixbssmn2(ti,initialTg4$alpha, initialTg4$beta, initialTg4$nu, initialTg4$pii, g=4, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g5    <- EMmixbssmn2(ti,initialTg5$alpha, initialTg5$beta, initialTg5$nu, initialTg5$pii, g=5, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_t_g6    <- EMmixbssmn2(ti,initialTg6$alpha, initialTg6$beta, initialTg6$nu, initialTg5$pii, g=6, "t", accuracy = 10^-6, iter.max = 500,aitken=TRUE)


round(c(fit_t_g1$result$lk   ,fit_t_g2$result$lk   ,fit_t_g3$result$lk   , fit_t_g4$result$lk),digits = 4)
c(fit_t_g1$result$aic  ,fit_t_g2$result$aic  ,fit_t_g3$result$aic  , fit_t_g4$result$aic)
c(fit_t_g1$result$bic  ,fit_t_g2$result$bic  ,fit_t_g3$result$bic  , fit_t_g4$result$bic)
c(fit_t_g1$result$iter,fit_t_g2$result$iter,fit_t_g3$result$iter, fit_t_g4$result$iter)
round(c(fit_t_g1$result$rates,fit_t_g2$result$rates,fit_t_g3$result$rates, fit_t_g4$result$rates),digits = 4)

round(InfmatrixmixBSSMN1(ti,fit_t_g2$result$pii,fit_t_g2$result$alpha,fit_t_g2$result$beta,fit_t_g2$result$nu,mfamily="t")$EP,digits = 4)

#-----------------------------------------------------------------------------------------------------
initialSLg1  <- initial.Values(ti,g=1,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialSLg2  <- initial.Values(ti,g=2,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = TRUE)
initialSLg3  <- initial.Values(ti,g=3,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialSLg4  <- initial.Values(ti,g=4,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialSLg5  <- initial.Values(ti,g=5,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialSLg6  <- initial.Values(ti,g=6,algorithm="k-means","sl",lower=1,upper=20,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)


fit_sl_g1  <- EMmixbssmn2(ti,initialSLg1$alpha, initialSLg1$beta, initialSLg1$nu, initialSLg1$pii, g=1, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE,rates=FALSE)
fit_sl_g2  <- EMmixbssmn2(ti,initialSLg2$alpha, initialSLg2$beta, initialSLg2$nu, initialSLg2$pii, g=2, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_sl_g3  <- EMmixbssmn2(ti,initialSLg3$alpha, initialSLg3$beta, initialSLg3$nu, initialSLg3$pii, g=3, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_sl_g4  <- EMmixbssmn2(ti,initialSLg4$alpha, initialSLg4$beta, initialSLg4$nu, initialSLg4$pii, g=4, "sl", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_sl_g5  <- EMmixbssmn2(ti,initialSLg5$alpha, initialSLg5$beta, initialSLg5$nu, initialSLg5$pii, g=5, "sl", accuracy = 10^-6, iter.max = 1000,aitken=TRUE)
fit_sl_g6  <- EMmixbssmn2(ti,initialSLg6$alpha, initialSLg6$beta, initialSLg6$nu, initialSLg6$pii, g=6, "sl", accuracy = 10^-6, iter.max = 1000,aitken=TRUE)


c(fit_sl_g1$result$lk,fit_sl_g2$result$lk,fit_sl_g3$result$lk,fit_sl_g4$result$lk)
c(fit_sl_g1$result$aic,fit_sl_g2$result$aic,fit_sl_g3$result$aic,fit_sl_g4$result$aic)
c(fit_sl_g1$result$bic,fit_sl_g2$result$bic,fit_sl_g3$result$bic,fit_sl_g4$result$bic)
c(fit_sl_g1$result$iter,fit_sl_g2$result$iter,fit_sl_g3$result$iter,fit_sl_g4$result$iter)
c(fit_sl_g1$result$rates,fit_sl_g2$result$rates,fit_sl_g3$result$rates,fit_sl_g4$result$rates)

round(InfmatrixmixBSSMN1(ti,fit_sl_g2$result$pii,fit_sl_g2$result$alpha,fit_sl_g2$result$beta,fit_sl_g2$result$nu,mfamily="t")$EP,digits = 4)


#-----------------------------------------------------------------------------------------------------
initialCNg1  <- initial.Values(ti,g=1,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialCNg2  <- initial.Values(ti,g=2,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = TRUE)
initialCNg3  <- initial.Values(ti,g=3,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialCNg4  <- initial.Values(ti,g=4,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialCNg5  <- initial.Values(ti,g=5,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
initialCNg6  <- initial.Values(ti,g=6,algorithm="k-means","cn",lower=0.1,upper=0.9,space=0.01,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)


fit_cn_g1  <- EMmixbssmn2(ti,initialCNg1$alpha, initialCNg1$beta, initialCNg1$nu, initialCNg1$pii, g=1, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE,rates=FALSE)
fit_cn_g2  <- EMmixbssmn2(ti,initialCNg2$alpha, initialCNg2$beta, initialCNg2$nu, initialCNg2$pii, g=2, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_cn_g3  <- EMmixbssmn2(ti,initialCNg3$alpha, initialCNg3$beta, initialCNg3$nu, initialCNg3$pii, g=3, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_cn_g4  <- EMmixbssmn2(ti,initialCNg4$alpha, initialCNg4$beta, initialCNg4$nu, initialCNg4$pii, g=4, "cn", accuracy = 10^-6, iter.max = 500,aitken=TRUE)
fit_cn_g5  <- EMmixbssmn2(ti,initialCNg5$alpha, initialCNg5$beta, initialCNg5$nu, initialCNg5$pii, g=5, "cn", accuracy = 10^-6, iter.max = 1000,aitken=TRUE)
fit_cn_g6  <- EMmixbssmn2(ti,initialCNg6$alpha, initialCNg6$beta, initialCNg6$nu, initialCNg6$pii, g=6, "cn", accuracy = 10^-6, iter.max = 1000,aitken=TRUE)


InfmatrixmixBSSMN1(ti,fit_cn_g2$result$pii,fit_cn_g2$result$alpha,fit_cn_g2$result$beta,fit_cn_g2$result$nu,mfamily="cn")

c(fit_cn_g1$result$lk,fit_cn_g2$result$lk,fit_cn_g3$result$lk,fit_cn_g4$result$lk)
c(fit_cn_g1$result$aic,fit_cn_g2$result$aic,fit_cn_g3$result$aic,fit_cn_g4$result$aic)
c(fit_cn_g1$result$bic,fit_cn_g2$result$bic,fit_cn_g3$result$bic,fit_cn_g4$result$bic)
c(fit_cn_g1$result$iter,fit_cn_g2$result$iter,fit_cn_g3$result$iter,fit_cn_g4$result$iter)
c(fit_cn_g1$result$rates,fit_cn_g2$result$rates,fit_cn_g3$result$rates,fit_sl_g4$result$rates)

#Tabla de Estimativas
Estimates <- round(c(fitBHg2$result$alpha, fitBHg2$result$beta, fitBHg2$result$pii[1]),digits=4)
EP        <- round(InfmatrixmixBS2(ti,fitBHg2$result$pii,fitBHg2$result$alpha,fitBHg2$result$beta)$EP,digits = 4)
U         <- round(Estimates + 1.96*EP,digits=4)
L         <- round(Estimates - 1.96*EP,digits=4)


Estimates <- round(c(fit_cn_g2$result$alpha, fit_cn_g2$result$beta, fit_cn_g2$result$pii[1]),digits=4)
EP        <- round(InfmatrixmixBSSMN1(ti,fit_cn_g2$result$pii,fit_cn_g2$result$alpha,fit_cn_g2$result$beta,fit_cn_g2$result$nu,mfamily="cn")$EP,digits=4)
U         <- round(Estimates + 1.96*EP,digits=4)
L         <- round(Estimates - 1.96*EP,digits=4)


Estimates <- round(c(fit_t_g2$result$alpha, fit_t_g2$result$beta, fit_t_g2$result$pii[2]),digits=4)
EP        <- round(InfmatrixmixBSSMN1(ti,fit_t_g2$result$pii,fit_t_g2$result$alpha,fit_t_g2$result$beta,fit_t_g2$result$nu,mfamily="t")$EP,digits=4)
U         <- round(Estimates + 1.96*EP,digits=4)
L         <- round(Estimates - 1.96*EP,digits=4)

Estimates <- round(c(fit_sl_g2$result$alpha, fit_sl_g2$result$beta, fit_sl_g2$result$pii[2]),digits=4)
EP        <- round(InfmatrixmixBSSMN1(ti,fit_sl_g2$result$pii,fit_sl_g2$result$alpha,fit_sl_g2$result$beta,fit_sl_g2$result$nu,mfamily="sl")$EP,digits=4)
U         <- round(Estimates + 1.96*EP,digits=4)
L         <- round(Estimates - 1.96*EP,digits=4)


postscript("MIX.FM.BS.SMN.densityBMI-NORMAL-g=1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="Normal")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg1$result$pii, fitBHg1$result$alpha, fitBHg1$result$beta), col="red", lty=1,lwd=2)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityBMI-T-g=1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="T de student")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bssmn(xx, fit_t_g1$result$pii, fit_t_g1$result$alpha, fit_t_g1$result$beta, fit_t_g1$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityBMI-CN-g=1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="Normal contaminada")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bssmn(xx, fit_cn_g1$result$pii, fit_cn_g1$result$alpha, fit_cn_g1$result$beta, fit_cn_g1$result$nu, mfamily = "cn"), col="black", lty=1,lwd=2)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityBMI-SL-g=1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="Slash")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bssmn(xx, fit_sl_g1$result$pii, fit_sl_g1$result$alpha, fit_sl_g1$result$beta, fit_sl_g1$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
dev.off() #Fechando o dispositivo potscript




postscript("MIX.FM.BS.SMN.densityBMI-g2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="g=2")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg2$result$pii, fitBHg2$result$alpha, fitBHg2$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g2$result$pii, fit_t_g2$result$alpha, fit_t_g2$result$beta, fit_t_g2$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g2$result$pii, fit_cn_g2$result$alpha, fit_cn_g2$result$beta, fit_cn_g2$result$nu, mfamily = "cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g2$result$pii, fit_sl_g2$result$alpha, fit_sl_g2$result$beta, fit_sl_g2$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityBMI-g3eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="g=3")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg3$result$pii, fitBHg3$result$alpha, fitBHg3$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g3$result$pii, fit_t_g3$result$alpha, fit_t_g3$result$beta, fit_t_g3$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g3$result$pii, fit_cn_g3$result$alpha, fit_cn_g3$result$beta, fit_cn_g3$result$nu, mfamily = "cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g3$result$pii, fit_sl_g3$result$alpha, fit_sl_g3$result$beta, fit_sl_g3$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityBMI-g4eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="g=4")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg4$result$pii, fitBHg4$result$alpha, fitBHg4$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g4$result$pii, fit_t_g4$result$alpha, fit_t_g4$result$beta, fit_t_g4$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g4$result$pii, fit_cn_g4$result$alpha, fit_cn_g4$result$beta, fit_cn_g4$result$nu, mfamily = "cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g4$result$pii, fit_sl_g4$result$alpha, fit_sl_g4$result$beta, fit_sl_g4$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("MIX.FM.BS.SMN.densityBMI-g5eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
hist(ti, freq = FALSE, breaks = 22, xlab="t",col = "grey", border=FALSE,ylim=c(0,0.1),main="g=5")
xx  <- (seq(min(ti)-(max(ti)-min(ti))/(length(ti)-1),max(ti)+9*(max(ti)-min(ti))/(length(ti)-1),(max(ti)-min(ti))/(length(ti)-1)))
lines(xx, d.mixed.bs(xx, fitBHg5$result$pii, fitBHg5$result$alpha, fitBHg5$result$beta), col="red", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_t_g5$result$pii, fit_t_g5$result$alpha, fit_t_g5$result$beta, fit_t_g5$result$nu, mfamily = "t"), col="blue", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_cn_g5$result$pii, fit_cn_g5$result$alpha, fit_cn_g5$result$beta, fit_cn_g5$result$nu, mfamily = "cn"), col="black", lty=1,lwd=2)
lines(xx, d.mixed.bssmn(xx, fit_sl_g5$result$pii, fit_sl_g5$result$alpha, fit_sl_g5$result$beta, fit_sl_g5$result$nu, mfamily = "sl"), col="green", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript



postscript("cdfsimmixBMI_g1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=1",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,70,0.1)
df1    <- p.mixed.bs(xx,pii=fitBHg1$result$pii,alpha=fitBHg1$result$alpha,beta=fitBHg1$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g1$result$pii,alpha=fit_t_g1$result$alpha,beta=fit_t_g1$result$beta,nu = fit_t_g1$result$nu,mfamily="t")
df3    <- p.mixed.bssmn(xx,pii=fit_cn_g1$result$pii,alpha=fit_cn_g1$result$alpha,beta=fit_cn_g1$result$beta,nu = fit_cn_g1$result$nu,mfamily="cn")
df4    <- p.mixed.bssmn(xx,pii=fit_sl_g1$result$pii,alpha=fit_sl_g1$result$alpha,beta=fit_sl_g1$result$beta,nu = fit_sl_g1$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
lines(xx, df3,col="black", lty=1,lwd=2)
lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("cdfsimmixBMI_g2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=2",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,70,0.1)
df1    <- p.mixed.bs(xx,pii=fitBHg2$result$pii,alpha=fitBHg2$result$alpha,beta=fitBHg2$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g2$result$pii,alpha=fit_t_g2$result$alpha,beta=fit_t_g2$result$beta,nu = fit_t_g2$result$nu,mfamily="t")
df3    <- p.mixed.bssmn(xx,pii=fit_cn_g2$result$pii,alpha=fit_cn_g2$result$alpha,beta=fit_cn_g2$result$beta,nu = fit_cn_g2$result$nu,mfamily="cn")
df4    <- p.mixed.bssmn(xx,pii=fit_sl_g2$result$pii,alpha=fit_sl_g2$result$alpha,beta=fit_sl_g2$result$beta,nu = fit_sl_g2$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
lines(xx, df3,col="black", lty=1,lwd=2)
lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("cdfsimmixBMI_g3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=3",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,70,0.1)
df1    <- p.mixed.bs(xx,pii=fitBHg3$result$pii,alpha=fitBHg3$result$alpha,beta=fitBHg3$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g3$result$pii,alpha=fit_t_g3$result$alpha,beta=fit_t_g3$result$beta,nu = fit_t_g2$result$nu,mfamily="t")
df3    <- p.mixed.bssmn(xx,pii=fit_cn_g3$result$pii,alpha=fit_cn_g3$result$alpha,beta=fit_cn_g3$result$beta,nu = fit_cn_g2$result$nu,mfamily="cn")
df4    <- p.mixed.bssmn(xx,pii=fit_sl_g3$result$pii,alpha=fit_sl_g3$result$alpha,beta=fit_sl_g3$result$beta,nu = fit_sl_g2$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
lines(xx, df3,col="black", lty=1,lwd=2)
lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("cdfsimmixBMI_g4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=4",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,70,0.1)
df1    <- p.mixed.bs(xx,pii=fitBHg4$result$pii,alpha=fitBHg4$result$alpha,beta=fitBHg4$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g4$result$pii,alpha=fit_t_g4$result$alpha,beta=fit_t_g4$result$beta,nu = fit_t_g4$result$nu,mfamily="t")
df3    <- p.mixed.bssmn(xx,pii=fit_cn_g4$result$pii,alpha=fit_cn_g4$result$alpha,beta=fit_cn_g4$result$beta,nu = fit_cn_g4$result$nu,mfamily="cn")
df4    <- p.mixed.bssmn(xx,pii=fit_sl_g4$result$pii,alpha=fit_sl_g4$result$alpha,beta=fit_sl_g4$result$beta,nu = fit_sl_g4$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
lines(xx, df3,col="black", lty=1,lwd=2)
lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("cdfsimmixBMI_g5.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot.ecdf(ti, main="g=5",xlab="t",ylab="cdf",col="gray", lty=6)
xx     <- seq(0.01,70,0.1)
df1    <- p.mixed.bs(xx,pii=fitBHg5$result$pii,alpha=fitBHg5$result$alpha,beta=fitBHg5$result$beta)
df2    <- p.mixed.bssmn(xx,pii=fit_t_g5$result$pii,alpha=fit_t_g5$result$alpha,beta=fit_t_g5$result$beta,nu = fit_t_g5$result$nu,mfamily="t")
df3    <- p.mixed.bssmn(xx,pii=fit_cn_g5$result$pii,alpha=fit_cn_g5$result$alpha,beta=fit_cn_g5$result$beta,nu = fit_cn_g5$result$nu,mfamily="cn")
df4    <- p.mixed.bssmn(xx,pii=fit_sl_g5$result$pii,alpha=fit_sl_g5$result$alpha,beta=fit_sl_g5$result$beta,nu = fit_sl_g5$result$nu,mfamily="sl")
lines(xx, df1,col="blue", lty=2,lwd=2)
lines(xx, df2,col="red", lty=3,lwd=2)
lines(xx, df3,col="black", lty=1,lwd=2)
lines(xx, df4,col="brown", lty=1,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript


library(survival)

fit     <- survfit(Surv(ti) ~ 1)
postscript("survival_apliBMIg1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=1")
xx     <- seq(0.1,length(ti),0.1)
#df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg1$result$pii,alpha=fitBHg1$result$alpha,beta=fitBHg1$result$beta)
#df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g1$result$pii,alpha=fit_t_g1$result$alpha,beta=fit_t_g1$result$beta,nu=fit_t_g1$result$nu,mfamily = "t")
#df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g1$result$pii,alpha=fit_cn_g1$result$alpha,beta=fit_cn_g1$result$beta,nu=fit_cn_g1$result$nu,mfamily = "cn")
#df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g1$result$pii,alpha=fit_sl_g1$result$alpha,beta=fit_sl_g1$result$beta,nu=fit_sl_g1$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
lines(xx, df_g_3,col="black", lty=1,lwd=2)
lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("survival_apliBMIg2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=2")
xx     <- seq(0.1,length(ti),0.1)
df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg2$result$pii,alpha=fitBHg2$result$alpha,beta=fitBHg2$result$beta)
df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g2$result$pii,alpha=fit_t_g2$result$alpha,beta=fit_t_g2$result$beta,nu=fit_t_g2$result$nu,mfamily = "t")
df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g2$result$pii,alpha=fit_cn_g2$result$alpha,beta=fit_cn_g2$result$beta,nu=fit_cn_g2$result$nu,mfamily = "cn")
df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g2$result$pii,alpha=fit_sl_g2$result$alpha,beta=fit_sl_g2$result$beta,nu=fit_sl_g2$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
lines(xx, df_g_3,col="black", lty=1,lwd=2)
lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("survival_apliBMIg3.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=3")
xx     <- seq(0.1,length(ti),0.1)
df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg3$result$pii,alpha=fitBHg3$result$alpha,beta=fitBHg3$result$beta)
df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g3$result$pii,alpha=fit_t_g3$result$alpha,beta=fit_t_g3$result$beta,nu=fit_t_g3$result$nu,mfamily = "t")
df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g3$result$pii,alpha=fit_cn_g3$result$alpha,beta=fit_cn_g3$result$beta,nu=fit_cn_g3$result$nu,mfamily = "cn")
df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g3$result$pii,alpha=fit_sl_g3$result$alpha,beta=fit_sl_g3$result$beta,nu=fit_sl_g3$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
lines(xx, df_g_3,col="black", lty=1,lwd=2)
lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript

postscript("survival_apliBMIg4.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
plot(fit, xlab="t",ylab="S(t)",main="g=4")
xx     <- seq(0.1,length(ti),0.1)
df_g_1 <- 1-p.mixed.bs(xx,pii=fitBHg4$result$pii,alpha=fitBHg4$result$alpha,beta=fitBHg4$result$beta)
df_g_2 <- 1-p.mixed.bssmn(xx,pii=fit_t_g4$result$pii,alpha=fit_t_g4$result$alpha,beta=fit_t_g4$result$beta,nu=fit_t_g4$result$nu,mfamily = "t")
df_g_3 <- 1-p.mixed.bssmn(xx,pii=fit_cn_g4$result$pii,alpha=fit_cn_g4$result$alpha,beta=fit_cn_g4$result$beta,nu=fit_cn_g4$result$nu,mfamily = "cn")
df_g_4 <- 1-p.mixed.bssmn(xx,pii=fit_sl_g4$result$pii,alpha=fit_sl_g4$result$alpha,beta=fit_sl_g4$result$beta,nu=fit_sl_g4$result$nu,mfamily = "sl")
lines(xx, df_g_1,col="blue",lty=2,lwd=2)
lines(xx, df_g_2,col="red", lty=2,lwd=2)
lines(xx, df_g_3,col="black", lty=1,lwd=2)
lines(xx, df_g_4,col="brown", lty=2,lwd=2)
legend("topright",c("N","T","NC","SL"),col=c("red", "blue","black", "green"),lty=c(2,2,1,2),lwd=2,bty="n",inset=0.1,cex=0.8)
dev.off() #Fechando o dispositivo potscript
