library(tidyverse)
library(foreach)
library(doParallel)
library(vegan)
library(VGAM)
library(MASS)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set Parameters, D=Dimensionality, M=Magnification Factor
D <- seq(5,200,by=50)
#M <- c(1:10,seq(12,98,by=2),seq(100,500,by=25),seq(550,1500,by=50),seq(2000,10000,by=500))
M <- c(1,2)

# Elaborates Transformations Effects
#------------------------------------------------------------------------------#
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl) 

set.seed(42)
df <-
    foreach(d=D, .combine="rbind") %:%
    foreach(m=M, .combine="rbind", .packages=c("ToyModel")) %dopar% {
      
      toy <- toy_model(n=10^4, cor=diag(d), M=m,
                       qdist=qnorm, 
                       param=c(mean=0, sd=1),
                       method="pearson",
                       force.positive=TRUE)
      
      data.frame("d"=d, "m"=m, 
                 "ERR_L1"= mean(abs(toy$cor_NorTA - toy$cor_L1)),
                 "ERR_CLR"=mean(abs(toy$cor_NorTA - toy$cor_CLR)),
                 "pielou"=mean(apply(toy$NorTA,1,vegan::diversity) / log(d)))
      }
stopCluster(cl)
saveRDS(df, "Compositional_Effects.rds")

