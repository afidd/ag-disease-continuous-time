library(survival)
library(KMsurv)
library(OIsurv)

x<-read.csv("susceptible.csv")
surv.object<-Surv(x$value, x$censored)
x.fit<-survfit(surv.object ~ 1)

pdf("cens_dist.pdf")
plot(x.fit)
dev.off()
