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

samples <- 50000 # samples to be kept after burn in
thin <- 10
N <- samples / thin

chain1 <- read.csv(file=paste("outputs/chain_iter_1_to_", samples,
                     "_chain1.csv", sep=""), header=TRUE)
chain2 <- read.csv(file=paste("outputs/chain_iter_1_to_", samples,
                     "_chain2.csv", sep=""), header=TRUE)
chain3 <- read.csv(file=paste("outputs/chain_iter_1_to_", samples,
                     "_chain3.csv", sep=""), header=TRUE)
chain4 <- read.csv(file=paste("outputs/chain_iter_1_to_", samples,
                     "_chain4.csv", sep=""), header=TRUE)

chain1 <- as.matrix(chain1)
chain2 <- as.matrix(chain2)
chain3 <- as.matrix(chain3)
chain4 <- as.matrix(chain4)


png(filename="plots/chains.png")

for (i in 1:6) {

  plot(c(1:N), chain1[,i+1], type='l', col='black',
       main=paste("alpha",i), xlab="iteration no.",
       ylab=paste("alpha",i))
  points(c(1:N), chain2[,i+1], type='l', col='blue')
  points(c(1:N), chain3[,i+1], type='l', col='green')
  points(c(1:N), chain4[,i+1], type='l', col='red')

}

dev.off()
