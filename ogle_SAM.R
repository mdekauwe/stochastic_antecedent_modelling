#!/usr/bin/Rscript

# Ogle et al.'s stochastic antecedent modelling (SAM) framework
#
# - Attempt to port Ogle's ANPP OpenBUGS example from Appendix 2. See Box 2 in
#   the main text.
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

wd <- getwd()
setwd(wd)

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

# creating the list of data to send to JAGS
data = list('block'=block, 'YearID'=YearID, 'Event'=Event, 'ppt'=ppt)

samples <- 10000
burn <- samples * 0.1
nadapt <- 100  # adaptions to tune sampler
nchains <- 4
# thinning rate, save every 10th iteration to reduce correlation between
# consecutive values in the chain
thin <- 10
jags <- jags.model('ogle_model.R', data=data, n.chains=nchains, n.adapt=nadapt)
mcmc_samples <- coda.samples(jags, n.iter=samples, n.burnin=burn, thin=thin,
                             variable.names=c('mu', "tau"))
#mcmc_samples <- coda.samples(jags, n.iter=samples, n.burnin=burn, thin=thin,
#                             variable.names=c('NPP', 'a', 'Event', 'mu', 'sig',
#                                              'tau', 'mu_ev', 'sig_ev', 'tau_ev'))

#plot(mcmc_samples)
summary(mcmc_samples)
