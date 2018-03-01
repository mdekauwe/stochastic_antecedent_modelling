#!/usr/bin/Rscript

# Plot the chains to check on mixing ...
#
# Author: Martin De Kauwe
# Email: mdekauwe@gmail.com
# Date: 01.03.2018

library(ggplot2)
library(cowplot)


wd <- getwd()
setwd(wd)

chain1 <- read.csv(file=paste("outputs/chain_iter_1_to_", chain,
                     "_chain1.csv", sep=""), header=TRUE)
chain2 <- read.csv(file=paste("outputs/chain_iter_1_to_", chain,
                     "_chain2.csv", sep=""), header=TRUE)
chain3 <- read.csv(file=paste("outputs/chain_iter_1_to_", chain,
                     "_chain3.csv", sep=""), header=TRUE)
chain4 <- read.csv(file=paste("outputs/chain_iter_1_to_", chain,
                     "_chain4.csv", sep=""), header=TRUE)

samples <- 50000 # samples to be kept after burn in

png(filename="plots/chains.png")

for (i in 1:6) {

  plot(c(1:samples), chain1[,i+1], type='l', col='black',
       main=paste("alpha",i), xlab="iteration no.",
       ylab=paste("alpha",i))
  points(c(1:samples), chain2[,i+1], type='l', col='blue')
  points(c(1:samples), chain3[,i+1], type='l', col='green')
  points(c(1:samples), chain4[,i+1], type='l', col='red')

}

dev.off()
