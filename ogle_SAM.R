#!/usr/bin/Rscript

# Ogle et al.'s stochastic antecedent modelling (SAM) framework
#
# Port of Ogle's ANPP OpenBUGS example from Appendix 2. See Box 2 in
# the main text.
#
# This is the wrapper function: sets things up, calls the model script and
# does something with the result
#
# Reference
# ---------
# * Ogle et al. (2015) Quantifying ecological memory in plant and ecosystem
#   processes. Ecology Letters, 18: 221–235
#
# Author: Martin De Kauwe
# Email: mdekauwe@gmail.com
# Date: 22.02.2018

library(rjags)
library(ggplot2)
library(cowplot)

wd <- getwd()
setwd(wd)

N <- 52
Nyrs <- 91
Nblocks <- 38

# Number of past years, including the current year for which the antecedent
# conditions are computed
Nlag <- 5
INCH_TO_MM <- 25.4

#
## Setup driving data ...
#

# dataset 1
#
# the time block that each month is assigned to such that for 60 different
# months, we are only estimating 38 unique monthly weights
block = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
          18, 19, 20, 21, 22, 23, 24, 25, 25, 26, 26, 27, 27, 28, 28,
          29, 29, 30, 30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34,
          34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38,
          38, 38)
block <- matrix(block, 5, 12)

# ANPP and precipitation event data for each year, extracted from Lauenroth
# and Sala (1992).
df2 = read.table("data/dataset2.csv", na.strings="NA", skip=1, sep=" ",
                 stringsAsFactors=FALSE, header=TRUE)
Event <- df2[c("Event1", "Event2", "Event3", "Event4")]
YearID <- df2$YearID

# Monthly precipitation data
# - data were obtained for Fort Collins, Colorado, which is located about
# 56 miles southwest of the CPER site. Fort Collins' monthly precipitation
# totals (mm) for 1893-2009 were downloaded from the Western Regional Climate
# Center. Ideally, we would want to use monthly precipitation from the study
# site (CPER), but this site did not provide complete monthly records for
# multiple years prior to the onset of the ANPP measurements.
df3 = read.table("data/dataset3.csv",  na.strings="NA", skip=1, sep=" ",
                 stringsAsFactors=FALSE, header=TRUE)
ppt <- df3[c("ppt1", "ppt2", "ppt3", "ppt4", "ppt5", "ppt6", "ppt7", "ppt8",
             "ppt9", "ppt10", "ppt11", "ppt12")]


#
## MCMC ...
#

# creating the list of data to send to JAGS
data = list('block'=block, 'YearID'=YearID, 'Event'=Event, 'ppt'=ppt,
            'Nlag'=Nlag, 'N'=N, 'Nyrs'=Nyrs, 'Nblocks'=Nblocks,
            'INCH_TO_MM'=INCH_TO_MM)

samples <- 100000 # samples to be kept after burn in
burn <- samples * 0.1 # iterations for burn in
nadapt <- 100  # adaptions to tune sampler
nchains <- 4
# thinning rate, save every 10th iteration to reduce correlation between
# consecutive values in the chain
thin <- 10
jags <- jags.model('model.R', data=data, n.chains=nchains, n.adapt=nadapt)
fit <- coda.samples(jags, n.iter=samples, n.burnin=burn, n.thin=thin,
                    variable.names=c('mu','alpha','deviance','Dsum'))

#
## Extract ouputs
#

pars <- summary(fit)$statistics[,1]

pred <- rep(NA, N)
offset <- 8
for (i in 1:N) {

  pred[i] <- pars[offset+i]

}

plot(df2$Year, pred, col="salmon", xlim=range(c(1940, 1990)),
     ylim=range(c(-700, 200)))
points(df2$Year, df2$NPP, col="royalblue")


#plot(fit, trace=FALSE, density=TRUE)
#plot(fit, trace=TRUE, density=FALSE)


for (i in 1:nchains) {

  # Extract fitted model
  write.csv(fit[[i]], file=paste("outputs/samples_iter_1_to_",
            samples, sprintf("_chain%i.csv", i), sep=""), row.names=FALSE)

  # Save states
  ss <- coef(jags, chain=i)
  save(ss, file=paste("outputs/saved_state_iter_1_to_", samples,
       sprintf("_chain%i.bin", i), sep=""))

}

# Before assessing the Gelman criteria, first check that hte posterior
# distributions are all approximately Normal.
densplot(fit)

# Rhat values for Gelman criteria
# NOT CURRENTLY WORKING - DON'T KNOW WHY.
g <- matrix(NA, nrow=nvar(fit), ncol=51)

for (v in 1:nvar(fit)) {

  # compare chains to check on mixing
  x <- gelman.plot(fit[,v])
  y <- x$shrink
  g[v,] <- y[,,1]

}

out <- rbind(x$last.iter - 100, g)
write.csv(t(out), file=paste("outputs/samples_iter_1_to_", samples,
          "_Rhat.csv", sep=""), row.names=FALSE)

# Once the MCMC has converged to the posterior distribution, we compute hte
# DIC by running the MCMC 1000 more iterations using:
DIC.calc <- dic.samples(jags, n.iter=1000, type="pD")
print(DIC.calc)
