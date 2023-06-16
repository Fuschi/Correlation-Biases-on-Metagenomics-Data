library(tidyverse)
library(ToyModel)
library(gridExtra)
library(SpiecEasi)


otu <- readRDS("data/otu_HMP2.rds")

# Filtering Rarest OTUs
#------------------------------------------------------------------------------#
prev <- function(X) colSums(X>0)/nrow(X)
medi <- function(X) apply(X,2,function(x)median(x[x>0]))

otu.filt <- otu[, prev(otu)>=.25]
otu.filt <- otu.filt[, medi(otu.filt)>=5]
otu.filt <- otu.filt[rowSums(otu.filt)>=10^3,]

# Set random seed for reproducibility
# However, if it is not set, the results remain very stable, only for 
# the maximum values of depth and counts will greater variations be observed
set.seed(42)

# Hurdle Truncated Log-Normal data generation
#------------------------------------------------------------------------------#
fit.hln <- t(apply(otu.filt,2,function(x)mle.htrlnorm(x)$estimate, simplify=TRUE))
fit.hln <- cbind(fit.hln,"b"=apply(otu.filt,2,function(x) log(max(x))))
count.hln <- NorTA(n=nrow(otu.filt), cor=cor(clr(otu.filt+1)), qdist=qhtrlnorm,
                   param=fit.hln)

# Zero inflated negative binomial
#------------------------------------------------------------------------------#
count.znb <- synth_comm_from_counts(otu.filt,dist="zinegbin")


summary.count <- rbind(
  quantile(otu.filt,probs=seq(0,1,by=.1)),
  quantile(count.znb,probs=seq(0,1,by=.1)),
  quantile(count.hln,probs=seq(0,1,by=.1))
)
rownames(summary.count) <- c("HMP2", "zinegbin","htrlnorm")
summary.count <- round(summary.count)

summary.depth <- rbind(
  quantile(rowSums(otu.filt),probs=seq(0,1,by=.1)),
  quantile(rowSums(count.znb),probs=seq(0,1,by=.1)),
  quantile(rowSums(count.hln),probs=seq(0,1,by=.1))
)
rownames(summary.depth) <- c("HMP2", "zinegbin","htrlnorm")
summary.depth <- round(summary.depth)