# This reads csv files of total outbreak sizes and makes
# a single plot showing survival curves of both together.
library(survival)
library(KMsurv)
library(OIsurv)

plot_total_outbreak <- function(inputa, inputb, outfile) {
  x<-read.csv(inputa)
  x$censor<-rep(1, times=length(x$outbreaksize))
  surv.object<-Surv(x$outbreaksize, x$censor)
  x.fit<-survfit(surv.object ~ 1)
  y<-read.csv(inputb)
  ysurv.object<-Surv(y$outbreaksize, rep(1, times=length(y$outbreaksize)))
  y.fit<-survfit(ysurv.object ~ 1)

  pdf(outfile)
  plot(y.fit)
  lines(x.fit, col="blue")
  title("Survival for Total Outbreak Size", xlab = "Count of Farms",
    ylab="Survival Fraction")
  dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
plot_total_outbreak(args[1], args[2], "total_outbreak_compare.pdf")
