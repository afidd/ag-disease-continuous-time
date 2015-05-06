# Makes a plot of the survival of a censored trial
# saved as a csv with two columns, "value" and "censored"
# where "censored" is 0 for censored, 1 for present.
# R --no-save --args "clinical" < individual_dist.R

library(survival)
library(KMsurv)
library(OIsurv)

plot_survival <- function(name) {
  x<-read.csv(paste(name, ".csv", sep=""))
  surv.object<-Surv(x$value, x$censored)
  x.fit<-survfit(surv.object ~ 1)

  pdf(paste(name, ".pdf", sep=""))
  plot(x.fit)
  title(paste("Survival for", name),
    xlab = "Days", ylab = "Survival Fraction")
  dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
plot_survival(args[1])
