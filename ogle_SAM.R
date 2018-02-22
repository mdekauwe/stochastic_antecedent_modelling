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
ppt = read.table("data/dataset3.csv",  na.strings="NA", skip=1, sep=" ",
                 stringsAsFactors=FALSE, header=TRUE)

#creating the list of data to send to JAGS
data <- list()
data[[1]] <- Trt
data[[2]] <- d2$GPP
data[[3]] <- flatten(SWC_DSM1)#
data[[4]] <- flatten(Temp_DST1)#
data[[5]] <- flatten(VPDh)
data[[6]] <- PARh
data[[7]] <- flatten(Gness)#
data[[8]] <- flatten(Abg_Nit)#
data[[9]] <- flatten(VPD_antX2)#
data[[10]] <- flatten(SWC_DSM1_antX2)
data[[11]] <- flatten(Temp_DST1_antX2)#
data[[12]] <- d2$year
data[[13]] <- d2$Plot_Ed
data[[14]] <- d2$Hour_GPP
data[[15]] <- d2$DayN
data[[16]] <- flatten(GnessDiff_antX2)   #
data[[17]] <- d2$CO2umol
data[[18]] <- CO2switch



n_adapt = 100
poodel = jags.model("ogle_model.R",
                    data=datalist,
                    inits=inits,
                    n.adapt=n_adapt,
                    n.chains=length(inits))
