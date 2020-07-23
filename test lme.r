# Broad and narrow sense inference and MA for linear model

library(nlme)
library(lme4)
library(metafor)
library(emmeans)

# example data used in Dixon and Moore book chapter (4 days)
#trigly <- read.table('triglybook.txt', header=T, as.is=T)

# trigly has 2 extra days with very different method effects

trigly <- read.table('trigly.txt', header=T, as.is=T)
trigly$method <- factor(trigly$method)
trigly$day <- factor(trigly$day)

triglyGrp <- groupedData(y ~ method | day, data=trigly)

trigly.lme1 <- lme(y ~ method, random = ~ 1, data=triglyGrp)
# | day implicit by the groupedData specification

trigly.lme1b <- lme(y ~ method, random = ~ 1|day, data=trigly)
# direct specification works here
summary(trigly.lme1)$tTable

# fit random interaction: only for group 2, 
#   but allow correl with intercept
trigly.lme2 <- update(trigly.lme1, 
  random = ~ method
  )
# note: update of trigly.lme1b doesn't work, 
#   need to start with grouped data specification

summary(trigly.lme2)$tTable

# see what's happening with the estimated coefficients
coef(trigly.lme1)
coef(trigly.lme2)

trigly.lmer1 <- lmer(y ~ method + (1|day), data=trigly)
trigly.lmer2 <- lmer(y ~ method + (method|day), data=trigly)

coef(trigly.lmer1)
coef(trigly.lmer2)
summary(trigly.lmer2)

# Repeat using meta analysis
# 2 ways to estimate day-specific difference and se

# fit model separately to each day
#   3 parameters per day (intercept, method, and error var.)
ndays <- length(unique(trigly$day))
trigly.ests <- matrix(NA, nrow=ndays, ncol=2)
for (i in 1:ndays) {
  bit <- trigly[trigly$day ==i,]
  temp <- summary(lm(y ~ method, data=bit))
  trigly.ests[i,] <- temp$coef[2,1:2]
}

# use lmList to fit all days
#   default is to use pooled error variance

temp <- lmList(y ~ method | day, data=triglyGrp)
trigly.ests2 <- summary(temp)$coef[,1:2,'method2']

# broad sense inference
rma(y=trigly.ests[,1], se=trigly.ests[,2], 
  method='REML')
rma(y=trigly.ests[,1], se=trigly.ests[,2], 
  method='REML', weighted=F)

# narrow sense inference
rma(y=trigly.ests[,1], se=trigly.ests[,2], 
  method='FE', weighted=F)

# narrow sense inference, done right with  fixed interaction
trigly.int <- lm(y ~ method*day, data=trigly)
trigly.emm <- emmeans(trigly.int, 'method')
pairs(trigly.emm)

# using pooled error variance
rma(y=trigly.ests2[,1], se=trigly.ests2[,2], 
  method='REML', weighted=F)
# exactly the same as lmer() results, 
#   weighing not relevant because all ests. have same se


# plots for talk
pdf(file='MA.pdf', height=6, width=8)
par(mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))
matplot(1:6, trigly.ests[,1]+outer(trigly.ests[,2], c(-2,0,2)),
  xlim=c(0.5,8), type='n',
  xlab='Day', ylab='Difference', xaxt='n')
points(1:6, trigly.ests[,1], pch=19, col=4)
segments(1:6, trigly.ests[,1]-2*trigly.ests[,2],
  1:6, trigly.ests[,1]+2*trigly.ests[,2] )
axis(1, at=1:6, labels=1:6)
abline(h=0, lty=3)
abline(v=6.5, lty=3)
points(c(7,7.5), rep(-8.09, 2), pch=18, cex=1.5, col=3)
segments(c(7,7.5), rep(-8.09, 2)-2*c(1.32, 5.94),
  c(7,7.5), rep(-8.09, 2)+2*c(1.32, 5.94) )
axis(1, at=c(7,7.5), labels=c('N','B'))
dev.off()
