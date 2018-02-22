#
# This is the ANPP model which is called by the wrapper script based on the
# example in Ogle et al.
#
# Reference
# ---------
# * Ogle et al. (2015) Quantifying ecological memory in plant and ecosystem
#   processes. Ecology Letters, 18: 221–235
#
# Author: Martin De Kauwe
# Email: mdekauwe@gmail.com
# Date: 22.02.2018

model {

  INCH_TO_MM <- 25.4
  N <- 52
  Nyrs <- 91
  Nblocks <- 38

  # Number of past years, including the current year for which the antecedent
  # conditions are computed
  Nlag <- 5

  #
  ## Compute anetcedent terms
  #

  # Dirichlet prior for monthly precipitation weights (due to restrictions
  # on when the built-in dirichlet distribution can be used, we are required
  # to use the relationship between the gamma distribution and the dirichlet
  # to assign the dirichlet prior. For each time block into the past, assign
  # the unnormalized weight (deltaX) a gamma(1,1) prior:
  #for (j in 1:Nblocks) {
  for (j in 1:38) {

    deltaX[j] ~ dgamma(1,1)

  }

  #for (t in 1:Nlag) {
  for (t in 1:5) {

    # Compute the yearly weights:
    yr_w[t] <- sum(weight[,t])
    alphad[t] <- 1

    for (m in 1:12) {

      # Redefine the unnormalized monthly weights to account for post-ANPP
      # harvest period; i.e., 2nd part involving equals and step functions
      # sets weight = 0 if in year 1 and in Oct, Nov, or Dec
      # (i.e., post- ANPP harvest).
      delta[m,t] <- (deltaX[block[t,m]]) * (1-equals(t,1) * step(m-9.5))

      # Compute normalized monthly weights, which will be between 0 and 1,
      # and will some to one.
      weight[m,t] <- delta[m,t] / sumD

      # Reorder the weights to go from most recent month (Dec of current
      # year) to “oldest” month (Jan at past year = Nlag).
      weightOrdered[(t-1)*12 + (12-m+1)] <- weight[m,t]

      # For each time into the past, compute the weighted precipitation
      # variable.
      #for (i in Nlag:Nyrs) {
      for (i in 5:91) {

        antX1[i,m,t] <- weight[m,t] * ppt[i-t+1,m]

      }

    }

  }

  # Compute sum of deltas (unnormalized weights), to be used to compute
  # the normalized antecedent weights:
  #for (t in 1:Nlag) {
  for (t in 1:5) {

    sumD1[t] <- sum(delta[,t])

  }
  sumD <- sum(sumD1[])

  # Compute the cumulative monthly weights:
  #for (t in 1:(12*Nlag)) {
  for (t in 1:(12*5)) {

    cum_weight[t] <- sum(weightOrdered[1:t])

  }

  # Compute the month within year weights (alpha’s = wP,m in Box 1 in main
  # text); that is, these weights sum to 1 within each past year
  for (m in 1:12) {

    #for (t in 1:Nlag) {
    for (t in 1:5) {

      alpha[m,t] <- delta[m,t] / sum(delta[,t])

    }

  }

  # Compute antecedent precipitation by summing the weighted precipitation
  # variable over months and past years ("exogenous variable"):
  #for (i in Nlag:Nyrs) {
  for (i in 5:91) {

    #for (t in 1:Nlag) {
    for (t in 1:5) {

      ant_sum1[i,t] <- sum(antX1[i,,t])

    }

    antX[i] <- sum(ant_sum1[i,])

  }

  #
  ## Model likelihood
  #

  #for (i in 1:N) {
  for (i in 1:52) {

    # Data model (or likelihood) for the observed NPP data:
    NPP[i] ~ dnorm(mu[i], tau)

    # Generate “replicated data” to evaluate model fit.
    NPP_rep[i] ~ dnorm(mu[i], tau)

    # Define model for latent (mean) NPP; Event[,k] represents the amount
    # of precipitation received in different size classes, where k indexes
    # the even size class (k=1 for < 5 mm; k=2 for 5-15 mm; k=3 for 15-
    # 30 mm; k=4 for >30 mm)
    mu[i] <- ( a[1] + (a[2] * antX[YearID[i]] * INCH_TO_MM) +
              (a[3] * Event[i,1]) + (a[4] * Event[i,2]) +
              (a[5] * Event[i,3]) + (a[6] * Event[i,4]) )

    # Compute first part of deviance
    D[i] <- log(2 * 3.1415926535) - log(tau) + (pow(NPP[i] - mu[i], 2) * tau)

    # Part of the calculation of the posterior predictive loss.
    # After each iteration we also compute Dsum, which is at end of this script.
    sq_diff[i] <- pow(NPP_rep[i] - NPP[i], 2)

    # Some of the precipitation event data are missing, so specify a simple
    # data model for the Event data for the purpose of estimating the
    # missing data:
    for (k in 1:4) {

      Event[i,k] ~ dnorm(mu_ev[k], tau_ev[k])

    }

  }

  #
  ## Priors
  #

  # Assign priors to the ANPP regression parameters. "a" is a vector of
  # coefficients that describes the effects of the exogenous (e.g. PPT) and
  # endogeneous (e.g. past events) covariates on mu.
  for (k in 1:6) {

    a[k] ~ dnorm(0, 1E-07)

  }

  # Prior for residual (observation) standard deviation, and compute
  # associated precision:
  sigma ~ dunif(0, 100)
  tau <- pow(sigma, -2) # precision (1 / variance)

  # Priors for parameters in the Event missing data model:
  for (k in 1:4) {

    mu_ev[k] ~ dunif(0, 500)
    sigma_ev[k] ~ dunif(0, 500)
    tau_ev[k] <- pow(sigma_ev[k], -2)

  }

  # Second part of deviance calculation
  deviance <- sum(D)

  # Posterior predictive loss is the posterior mean of Dsum
  Dsum <- sum(sq_diff[])

}
