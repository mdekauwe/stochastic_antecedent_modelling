#!/usr/bin/Rscript

# Plot the MCMC chains
#
# Author: Martin De Kauwe
# Email: mdekauwe@gmail.com
# Date: 22.02.2018

library(rjags)
library(ggplot2)
library(cowplot)

wd <- getwd()
setwd(wd)

nchains <- 4
samples <- 10000

samples1 <- read.csv(file=paste("outputs/samples_iter_1_to_", samples,
                     "_chain1.csv", sep=""), header=TRUE)
samples2 <- read.csv(file=paste("outputs/samples_iter_1_to_", samples,
                     "_chain2.csv", sep=""), header=TRUE)
samples3 <- read.csv(file=paste("outputs/samples_iter_1_to_", samples,
                     "_chain3.csv", sep=""), header=TRUE)
samples4 <- read.csv(file=paste("outputs/samples_iter_1_to_", samples,
                     "_chain4.csv", sep=""), header=TRUE)


end <- dim(samples1)[1]
start <- (end / 2) + 1
col = 8

# Calculate the deviance of the final 5000 iterations of each chain.
# If deviance diff > 5, then need to reinitialize.
mu_dev <- c(mean(samples1[start:end,col]), mean(samples2[start:end,col]),
            mean(samples3[start:end,col]), mean(samples4[start:end,col]))
min_dev <- c(min(samples1[start:end,col]), min(samples2[start:end,col]),
             min(samples3[start:end,col]), min(samples4[start:end,col]))
max_dev <- c(max(samples1[start:end,col]), max(samples2[start:end,col]),
            max(samples3[start:end,col]), max(samples4[start:end,col]))

print(mu_dev)
print(min_dev)
print(max_dev)

print(samples1["mu.1."])

hist(samples1["mu.1."])
#for (v in 1:nvar(samples1)) {
#
#  print(samples1[,v])
#
#
#}
