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

INCH_TO_MM <- 25.4
N <- 52
Nyrs <- 91
Nblocks <- 38

# Number of past years, including the current year for which the antecedent
# conditions are computed
Nlag <- 5

# the time block that each month is assigned to such that for 60 different
# months, we are only estimating 38 unique monthly weights
block = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
          18, 19, 20, 21, 22, 23, 24, 25, 25, 26, 26, 27, 27, 28, 28,
          29, 29, 30, 30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34,
          34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38,
          38, 38)
block <- matrix(block, 5, 12)


wd <- getwd()
setwd(wd)

# ANPP and precipitation event data for each year, extracted from Lauenroth
# and Sala (1992).
df2 = read.table("data/dataset2.csv", na.strings="NA", skip=1, sep=" ",
                 stringsAsFactors=FALSE, header=TRUE)

# Monthly precipitation data
# - data were obtained for Fort Collins, Colorado, which is located about
# 56 miles southwest of the CPER site. Fort Collins' monthly precipitation
# totals (mm) for 1893-2009 were downloaded from the Western Regional Climate
# Center. Ideally, we would want to use monthly precipitation from the study
# site (CPER), but this site did not provide complete monthly records for
# multiple years prior to the onset of the ANPP measurements.
df3 = read.table("data/dataset3.csv",  na.strings="NA", skip=1, sep=" ",
                 stringsAsFactors=FALSE, header=TRUE)

model {

  #
  ## Model likelihood
  #

  for (i in 1:N) {
    # Data model (or likelihood) for the observed NPP data:
    NPP[i] <- dnorm(mu[i], tau)

    # Generate “replicated data” to evaluate model fit.
    NPP_rep[i] <- dnorm(mu[i], tau)

    # Define model for latent (mean) NPP; Event[,k] represents the amount
    # of precipitation received in different size classes, where k indexes
    # the even size class (k=1 for < 5 mm; k=2 for 5-15 mm; k=3 for 15-
    # 30 mm; k=4 for >30 mm); convert antecedent precipitation (antX) from
    # inches to mm.
    mu[i] <- ( a[1] + (a[2] * antX[df2$YearID[i]] * INCH_TO_MM) +
              (a[3] * df2$Event[i,1]) + (a[4] * df2$Event[i,2]) +
              (a[5] * df2$Event[i,3]) + (a[6] * df2$Event[i,4]) )
  }



  # Compute antecedent precipitation by summing the weighted precipitation
  # variable over months and past years:
  for (i in Nlag:Nyrs) {
    for (t in 1:Nlag) {
      ant_sum1[i,t] <- sum(antX1[i,,t])
    }
    antX[i] <- sum(ant_sum1[i,])
  }


  #
  ## Priors
  #

  # Assign priors to the ANPP regression parameters (covariate effects):
  for (k in 1:6) {
    a[k] <- dnorm(0, 1E-07)
  }

  # Prior for residual (observation) standard deviation, and compute
  # associated precision:
  sig <- dunif(0, 100)
  tau <- pow(sig,-2)

  # Priors for parameters in the Event missing data model:
  for (k in 1:4) {
    mu_ev[k] <- dunif(0, 500)
    sig_ev[k] <- dunif(0, 500)
    tau_ev[k] <- pow(sig_ev[k], -2)
  }


}
