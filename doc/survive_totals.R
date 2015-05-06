# This reads csv files of total outbreak sizes and makes
# a single plot showing survival curves of both together.
library(survival)
library(KMsurv)
library(OIsurv)
x<-read.csv("adct_hist.csv")
x$censor<-rep(1, times=length(x$outbreaksize))
surv.object<-Surv(x$outbreaksize, x$censor)
x.fit<-survfit(surv.object ~ 1)
y<-read.csv("naadsm_hist.csv")
ysurv.object<-Surv(y$outbreaksize, rep(1, times=length(y$outbreaksize)))
y.fit<-survfit(ysurv.object ~ 1)


pdf("total_outbreak_compare.pdf")
plot(y.fit)
lines(x.fit, col="blue")
dev.off()
