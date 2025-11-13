library(tidyverse)
library(ToyModel)
library(VGAM)

# Populosity
n <- 10^3

# Correlation structure
D <- 5
R <- diag(D)
R[4,2] <- R[2,4] <- -.9

# Generation of Normal Multivariate distribution with desired correlation structure
rnorm <- mvtnorm::rmvnorm(n=n,sigma=R)

# Normal Quantile
unif <- stats::pnorm(rnorm)

# Transform of Marginal distribution
NorTA <- matrix(0,nrow=n, ncol=D)

NorTA[,1] <- qzinegbin(p=unif[,1], munb=20, size=30, pstr0=.25)
NorTA[,2] <- qzinegbin(p=unif[,2], munb=20, size=30, pstr0=.25)
NorTA[,3] <- qzinegbin(p=unif[,3], munb=20, size=30, pstr0=.25)
NorTA[,4] <- qzinegbin(p=unif[,4], munb=20, size=30, pstr0=.25)
NorTA[,5] <- qzinegbin(p=unif[,5], munb=20, size=30, pstr0=.25)

cor(NorTA[,2], NorTA[,4])
