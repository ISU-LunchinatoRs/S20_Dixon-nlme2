# gina.r: Philip's code for Gina's N leaching data set

library(dplyr)

library(nlme)
library(lme4)
library(nlraa)

library(emmeans)
library(metafor)

N <- read.csv('pro_lunchinators.csv', as.is=T)
N$rotation <- factor(N$rotation)

N$group <- with(N, paste(year, site_id, rotation, sep=':'))
N$siterot <- with(N, paste(year, site_id, sep=':'))

Ngrp <- groupedData(leaching_kgha ~ nrate_kgha | group, 
  outer = ~rotation,  
  data=N)

N10 <- subset(N, year==2010)
N10$siterot <- substring(N10$group, 6, 12)

Nklad <- subset(N, site_id=='klad')
Nklad$yearrot <- paste(Nklad$year, Nklad$rotation)

N10grp <- groupedData(leaching_kgha ~ nrate_kgha | siterot, 
  outer = ~rotation,  data=N10)
Nkladgrp <- groupedData(leaching_kgha ~ nrate_kgha | yearrot, 
#  outer = ~rotation, 
  data=Nklad)

pdf(file='N10.pdf', height=6, width=8)
plot(N10grp)
dev.off()


# fit bilinear curves to all groups
N.group <- nlsList(leaching_kgha ~ SSblin(nrate_kgha, Y0, b1, Yb, b2), 
  data=Ngrp)
  
# what do the subject-specific coefficients look like?
N.coef <- coefficients(N.group)
pairs(N.coef)

# 358 groups, 41 could not fit to single curve

# return to smaller data sets to demonstrate
N10.group <- nlsList(leaching_kgha ~ SSblin(nrate_kgha, Y0, b1, Yb, b2), 
  data=N10grp)
Nklad.group <- nlsList(leaching_kgha ~ SSblin(nrate_kgha, Y0, b1, Yb, b2), 
  data=Nkladgrp)

# Fit a mixed model to pool information across runs
N10.nlme <- nlme(N10.group)
Nklad.nlme <- nlme(Nklad.group)

plot(N10.nlme)
plot(augPred(N10.nlme))

plot(Nklad.nlme)
plot(augPred(Nklad.nlme))

nlsCoef <- coefficients(Nklad.group)
nlmeCoef <- coefficients(Nklad.nlme)

par(mar=c(3,3,0,0)+0.2, mgp=c(2, 0.8, 0))
par(mfrow=c(2,2))
for (i in 1:4) {
  plot(nlsCoef[,i], nlmeCoef[,i], pch=19, col=4,
    xlab='NLS estimate', ylab='NLME estimate')
  legend('topleft', bty='n', legend=names(nlsCoef)[i])
}

# Approach #1: treat rotation treatments as independent groups
# need to define when group the data

Ngrp <- groupedData(leaching_kgha ~ nrate_kgha | group, 
  outer= ~rotation,
  data=N)

N.nlme0 <- nlme(leaching_kgha ~ 
    SSblin(nrate_kgha, Y0, b1, Yb, b2), 
  fixed = list(Y0 + b1 + Yb + b2 ~ 1),
  random = Y0 + b1 + Yb + b2 ~ 1,
  data = Ngrp)
N.nlme0

# This should work, with either version, right?
N.nlme1 <- nlme(leaching_kgha ~ 
    SSblin(nrate_kgha, Y0, b1, Yb, b2), 
#  fixed = list(Y0 + b1 + Yb + b2 ~ rotation),
  fixed = list(Y0 ~ rotation, b1 ~ rotation, 
    Yb ~ rotation, b2 ~ rotation),
  random = Y0 + b1 + Yb + b2 ~ 1,
  data = Ngrp)
N.nlme1
# but result ignores rotation

# but this gets a bit further
N.nlme1 <- update(N.nlme0, 
  fixed = list(Y0 + b1 + Yb + b2 ~ rotation) )
# now looking for rotation-effect starting values
# can't get them from the self-start function

# have to remember how R codes factor effects
# use fixed effects from nlme0 as rotation 1 effect

N.fixef <- fixef(N.nlme0)

N.nlme1 <- update(N.nlme0, 
  fixed = list(Y0 + b1 + Yb + b2 ~ rotation),
  start = c(N.fixef[1], 0, N.fixef[2], 0,
    N.fixef[3],0,N.fixef[4],0)
   )

# if you had a lot of groups, here's a quick way to 
#   generate the vector of starting values 
c( rbind(N.fixef, rep(0, length(N.fixef))))

# if you had 5 levels of the factor, would have 4 copies
#   of the rep(), to generate 4 sets of 0's
# but only works when you want same # 0's for each parameter

summary(N.nlme1)
# looks like T statistic for Y0 rotation effect is small


N.nlme2 <- update(N.nlme0,
  fixed = list(Y0 ~ 1, b1 ~ rotation, 
    Yb ~ rotation, b2 ~ rotation),
  start = c(N.fixef[1], N.fixef[2], 0,
    N.fixef[3],0,N.fixef[4],0)
   )
# note only 7 starting values (only 1 value for Y0)

# to see differences between cc and cs, look at
#   rotationcs coefficients
summary(N.nlme2)
intervals(N.nlme2)
# note df: these are based on individuals, not groups

library(emmeans)
N.emm <- emmeans(N.nlme2, 'rotation', param='b2')
contrast(N.emm, 'pairwise')
# difference between emmeans and summary output may just be
#   round off error
# but still using # obs to get df
# emmeans is often "aware" of random effects
#   and uses kenward-rogers adjustments
# can do this for lme4 objects
#   but not lme objects

anova(N.nlme0, N.nlme2, N.nlme1)
# likelihood ratio test - legitimate here because models fit 
#   by ML (default for nlme)
# Can't compare fixed effect models when fit by REML

# An idea I'm toying with:
#   compare random effect sd's with and without covariates
#   tell you how much including covariate reduced variability

N.nlme0
N.nlme1

# Y0 almost no reduction; others more, especially Yb and b2

# Better approach: rotations are paired within a site/year combination

# requires multiple grouping levels 
#   site:year to indicate pairs, 
#   site:year:rotation to indicate each curve


# let's take a meta regression approach
# extract coefficients and standard errors

# first - repeat independent site/year/rotation analysis

temp <- summary(N.group)$coef
dim(temp)
names(temp[1,1,])
names(temp[1,,1])

# first column is the estimate, second is the se
#   extract those for each parameter
Y0 <- temp[,1:2, 'Y0']
b1 <- temp[,1:2, 'b1']
Yb <- temp[,1:2, 'Yb']
b2 <- temp[,1:2, 'b2']


# now have to figure out labels for the 305 entries
#  most reliable way is dplyr

temp <- N %>% group_by(group)
GrpLabel <- temp %>% summarize(N=n())
# GrpLabel$group is the label for each
#   doesn't matter what summary function you apply

GrpLabel <- GrpLabel %>% mutate(
  yearsite = substring(group, 1, 9), 
  rotation = substring(group, 11:12)
  )
# don't know why cs became s

GrpLabel$yearsite <- factor(GrpLabel$yearsite)
GrpLabel$rotation <- factor(GrpLabel$rotation)

# let's look at estimate vs se for each
# marking each by rotation

temp <- summary(N.group)$coef
par(mfrow=c(2,2), mar=c(3,3,0,0)+0.3, mgp=c(2,0.8,0))
for (i in 1:4) {
  plot(temp[,1,i], temp[,2,i], pch=19, 
    col=3+(GrpLabel$rotation=='cc'),
    xlab='Estimate', ylab='se')
  legend('topleft', bty='n', legend=names(temp[1,,1])[i])
}

for (i in 1:4) {
  boxplot(split(temp[,1,i], GrpLabel$rotation))
  legend('topleft', bty='n', legend=names(temp[1,1,])[i])
}

# focus on b2 - other parameters will use similar code

# fixed effects MA
rma(yi=b2[,1], sei=b2[,2], mods = ~rotation, 
  data=GrpLabel, method='FE', weighted=F)

# random effects meta regression for independent rotations model
rma(yi=b2[,1], sei=b2[,2], mods = ~rotation, 
  data=GrpLabel, method='REML', weighted =  F)

# account for blocking by siteYear
rma(yi=b2[,1], sei=b2[,2], mods = ~rotation + yearsite, 
  data=GrpLabel, method='REML')

# You could also construct difference between cc and cs "by hand"
#   the hard part will be getting the se
#   easy if assume independence, 
#   hard if account for correl within year/site

# what if you ignore the 'internal' se?
temp <- summary(N.group)$coef
head(temp[,1,])

# add all response variables to the data frame
GrpLabel <- cbind(GrpLabel, temp[,1,])

N.lmer <- lmer(b2 ~ rotation + (1 | yearsite), 
  data=GrpLabel) 
# this is narrow sense inference

# so, which "type" of MA is nlme fitting?
head(coef(N.nlme2))
# narrow sense inference = fixed effect MA!!

# Could fit the NLME equivalent of the random effect model 
#   i.e., random rotation:plot interaction 

N.nlme2 <- update(N.nlme0,
  random = pdDiag(Y0 + b1 + Yb + b2 ~ rotation) 
  )

# or

N.nlme1b <- update(N.nlme1,
  random = pdDiag(Y0 + b1 + Yb + b2 ~ rotation - 1) 
  )

# but 1st one takes > 30 minutes.  I halted it before it finished.
# I didn't try the second.

# here's trying with just the data from the klad site
# narrow sense inference
Nklad.nlme1 <- nlme(leaching_kgha ~ 
    SSblin(nrate_kgha, Y0, b1, Yb, b2), 
  fixed = list(Y0 + b1 + Yb + b2 ~ rotation),
  random = Y0 + b1 + Yb + b2 ~ 1,
  data = Nkladgrp)
# convergence problems.

# but this generates an error - can't figure out why
#   update() generates same error
Nklad.nlme2 <- nlme(leaching_kgha ~ 
    SSblin(nrate_kgha, Y0, b1, Yb, b2), 
  fixed = Y0 + b1 + Yb + b2 ~ rotation,
   random = Y0 + b1 + Yb + b2 ~ rotation, 
  data = Nkladgrp
  )
