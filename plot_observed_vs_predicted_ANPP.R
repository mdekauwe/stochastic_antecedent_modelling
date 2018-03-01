#!/usr/bin/Rscript

# Plot the predicted ANPP vs. the observed ANPP
#
# Author: Martin De Kauwe
# Email: mdekauwe@gmail.com
# Date: 01.03.2018

library(ggplot2)

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

en <- end(chain1)[1]
st <- en - N

# To use the MCMC samples for prediction, we combine the 4 chains into one.
mu_post <- rbind(chain1[st:en,9:60], chain2[st:en,9:60],
                 chain3[st:en,9:60], chain4[st:en,9:60])

mu_post_mean <- apply(mu_post, 2, mean)
mu_post_95CI <- apply(mu_post, 2, quantile, probs=c(0.025, 0.975))

lower <- mu_post_mean - mu_post_95CI[1,]
upper <-mu_post_95CI[2,] - mu_post_mean

# ANPP and precipitation event data for each year, extracted from Lauenroth
# and Sala (1992).
df2 = read.table("data/dataset2.csv", na.strings="NA", skip=1, sep=" ",
                 stringsAsFactors=FALSE, header=TRUE)
YearID <- df2$YearID
NPP <- df2$NPP


df <- data.frame(df2$Year, df2$NPP, mu_post_mean)
colnames(df) <- c("Year","NPP_obs", "NPP_pred")

ggplot(data=df, aes(x=Year, y=NPP_pred)) +
  geom_point(color="black") +
  geom_errorbar(aes(ymin=NPP_pred-lower, ymax=NPP_pred+upper),
                color="black") +
  geom_point(aes(x=Year, y=NPP_obs), colour="salmon") +
  xlab("Year") +
  ylab("NPP") +
  theme(aspect.ratio=1/1.618) +
  theme_classic()

ggsave("plots/observed_ANPP_vs_predicted.png", width=6, height=4)
