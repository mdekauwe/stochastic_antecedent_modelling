#!/usr/bin/env python

"""
Attempt to port Ogle's ANPP OpenBUGS example from Appendix 2. See Box 2 in the
main text.

Reference
---------
* Ogle et al. (2015) Quantifying ecological memory in plant and ecosystem
  processes. Ecology Letters, 18: 221–235
"""

import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import theano.tensor as tt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (23.12.2017)"
__email__   = "mdekauwe@gmail.com"

N = 52
# Number of past years, including the current year for which the antecedent
# conditions are computed
Nlag = 5
Nyrs = 91
Nblocks = 38

# the time block that each month is assigned to such that for 60 different
# months, we are only estimating 38 unique monthly weights
block = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,\
                  18, 19, 20, 21, 22, 23, 24, 25, 25, 26, 26, 27, 27, 28, 28,\
                  29, 29, 30, 30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34,\
                  34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38,\
                  38, 38]).reshape(5,12)


# ANPP and precipitation event data for each year, extracted from Lauenroth
# and Sala (1992).
df2 = pd.read_csv("data/dataset2.csv", na_values="NA", skiprows=1, sep=" ")

# Monthly precipitation data; data were obtained for Fort Collins, Colorado,
# which is located about 56 miles southwest of the CPER site. Fort Collins'
# monthly precipitation totals (mm) for 1893-2009 were downloaded from the
# Western Regional Climate Center. Ideally, we would want to use monthly
# precipitation from the study site (CPER), but this site did not provide
# complete monthly records for multiple years prior to the onset of the
# ANPP measurements.
df3 = pd.read_csv("data/dataset3.csv", na_values="NA", skiprows=1, sep=" ")

# https://github.com/takluyver/pymc/blob/master/pymc/sandbox/parse_winbugs.py
step = lambda x: x>0

delta = np.zeros((12,Nlag))

with pm.Model() as model:
    # Assign priors to the ANPP regression parameters (covariate effects)
    a = pm.Normal('a', mu=0, sd=1E-07, shape=6)

    # Prior for residual (observation) standard deviation, and compute
    # associated precision
    sig = pm.Uniform('sig', 0, 100)
    tau = tt.pow(sig, -2)

    # Priors for parameters in the Event missing data model:
    mu_ev = pm.Uniform('mu_ev', 0, 500, shape=4)
    sig_ev = pm.Uniform('sig_ev', 0, 500, shape=4)
    tau_ev = tt.pow(sig_ev, -2)

    # Some of the precipitation event data are missing, so specify a simple
    # data model for the Event data for the purpose of estimating the
    # missing data:
    Event = pm.Normal('Event', mu=mu_ev, tau=tau_ev, shape=4)

    # Dirichlet prior for monthly precipitation weights (due to restrictions
    # on when the built-in dirichlet distribution can be used, we are required
    # to use the relationship between the gamma distribution and the dirichlet
    # to assign the dirichlet prior. For each time block into the past, assign
    # the unnormalized weight (deltaX) a gamma(1,1) prior:
    deltaX = pm.Gamma('deltaX', 1, 1, shape=Nblocks)

    for t in range(Nlag):
        for m in range(12):
            # Redefine the unnormalized monthly weights to account for
            # post-ANPP # harvest period; i.e., 2nd part involving equals and
            # step functions # sets weight = 0 if in year 1 and in Oct, Nov,
            # or Dec (i.e., post- # ANPP harvest).
            delta[m,t] = (deltaX[block[t,m]]) * \
                            (1 - np.equal(t,1) * step(m - 9.5))

            # Compute normalized monthly weights, which will be between
            # 0 and 1, and will some to one.
            weight[m,t] = delta[m,t] / sumD

            # Reorder the weights to go from most recent month (Dec of current
            # year) to “oldest” month (Jan at past year = Nlag).
            weightOrdered[t*12 + (12-m+1)] = weight[m,t]

            # For each time into the past, compute the weighted precipitation
            # variable.
            for i in range(Nlag, Nyrs+1):
                antX1[i,m,t] = weight[m,t] * df3.ppt[i-t+1,m]

        # Compute the yearly weights:
        yr_w[t] = sum(weight[:,t])
        alphad[t] = 1

    # Compute sum of deltas (unnormalized weights), to be used to compute
    # the normalized antecedent weights:
    for t in range(Nlag):
        sumD1[t] = sum(delta[:,t])
    sumD = sum(sumD1)

    # Compute the cumulative monthly weights:
    for t in range(Nlag*12):
        cum_weight[t] = sum(weightOrdered[0:t])

    # Compute the month within year weights (alpha’s = wP,m in Box 1 in main
    # text); that is, these weights sum to 1 within each past year
    for m in range(12):
        for t in range(Nlag):
            alpha[m,t] = delta[m,t] / sum(delta[:,t])

    # Compute antecedent precipitation by summing the weighted precipitation
    # variable over months and past years:
    for i in range(Nlag, Nyrs+1):
        for t in range(Nlag):
            ant_sum1[i,t] = sum(antX1[i,:,t])
        antX[i] = sum(ant_sum1[i,:])

    for i in range(N):
        # Define model for latent (mean) NPP; Event[,k] represents the amount
        # of precipitation received in different size classes, where k indexes
        # the even size class (k=1 for < 5 mm; k=2 for 5-15 mm; k=3 for 15-
        # 30 mm; k=4 for >30 mm); convert antecedent precipitation (antX) from
        # inches to mm.
        mu[i] = a[0] + (a[1] * antX[df2.YearID[i]] * 25.4) + \
                (a[2] * df2.Event[i,0]) + (a[3] * df2.Event[i,1]) + \
                (a[4] * df2.Event[i,2]) + (a[5] * df2.Event[i,3])

        # Data model (or likelihood) for the observed NPP data:
        NPP[i] = pm.Normal('NPP', mu=mu[i], tau=tau)

        # Generate “replicated data” to evaluate model fit.
        NPP_rep[i] = pm.Normal('NPP_rep', mu=mu[i], tau=tau)


pm.traceplot(traces)
