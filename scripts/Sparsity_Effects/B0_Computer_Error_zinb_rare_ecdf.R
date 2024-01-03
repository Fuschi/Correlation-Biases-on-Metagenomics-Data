library(tidyverse)
library(doSNOW)
library(progress)
library(foreach)
library(VGAM)
library(SpiecEasi)
library(ToyModel)

setwd("~/Scrivania/Correlation-Biases-on-Metagenomics-Data/scripts/Sparsity_Effects/")

# Read HMP2library(tidyverse)
library(doSNOW)
otu <- readRDS("../../data/otu_HMP2.rds")
meta <- readRDS("../../data/meta_HMP2.rds")

# Select samples belonging to 69-001 subject in health status
otu.69001.H <- otu[meta$SubjectID=="69-001" & meta$CL4_2=="Healthy", ]

# Remove rarest OTUs using prevalence and median of non-zero values
otu.filt <- otu.69001.H[, colSums(otu.69001.H>0)/nrow(otu.69001.H) >= .33]
otu.filt <- otu.filt[, apply(otu.filt,2,function(x) median(x[x>0])>=5)]
otu.filt <- round(otu.filt/rowSums(otu.filt) * min(rowSums(otu.filt)))
rm(otu, otu.69001.H, meta)

# Fit ZINB params from real data for each filtered OTUs and
# save the first and ninen-th decil of the distribution
HMP2.params <- apply(otu.filt,2,function(x){
  SpiecEasi::fitdistr(as.numeric(x),"zinegbin")$par
}) %>% t %>% as.data.frame
HMP2.quantile.params <- HMP2.params %>% apply(2,function(x) quantile(x,probs=c(.1,.9)))
HMP2.params.filt <- HMP2.params %>%
  filter(munb>=HMP2.quantile.params["10%","munb"],
         munb<=HMP2.quantile.params["90%","munb"],
         size>=HMP2.quantile.params["10%","size"],
         size<=HMP2.quantile.params["90%","size"],
         pstr0>=HMP2.quantile.params["10%","pstr0"],
         pstr0<=HMP2.quantile.params["90%","pstr0"])


# Create the progression bar
nIteration <- 100
pb <- progress_bar$new(
  format = "[:bar] :elapsed | eta: :eta",
  total = nIteration*380,
  width = 60)
progress <- function(n){pb$tick()}

set.seed(42)
result <- data.frame()
for(iter in 1:nIteration){
  # START OUTER LOOP
  #----------------------------------------------------------------------------#
  
  # GENERATION OF THE UNDERLINE DATASET
  params_random_HMP2 <- data.frame(
    "munb"=quantile(HMP2.params.filt$munb, probs=runif(runif(200))),
    "size"=quantile(HMP2.params.filt$size, probs=runif(runif(200))),
    "pstr0"=quantile(HMP2.params.filt$pstr0, probs=runif(runif(200))))
  
  random_HMP2 <- ToyModel::toy_model(n=10^4, cor=diag(200), M=1,
                                     qdist=VGAM::qzinegbin,
                                     param=params_random_HMP2)
  random_cor0_HMP2 <- random_HMP2$cor_normal
  
  # CREATE PARAMETERS OF THE SUPERVISED COUPLE OF VARIABLES
  params_set <- expand_grid(
    "cor"=seq(-.9,.9,by=.1),
    "pstr0_1"=seq(0,.95,by=.05)
  ) %>%
    mutate(pstr0_2=pstr0_1, .after=pstr0_1) %>%
    mutate("munb_1"=quantile(HMP2.params.filt$munb, probs=runif(n())),
           "munb_2"=quantile(HMP2.params.filt$munb, probs=runif(n())),
           "size_1"=quantile(HMP2.params.filt$size, probs=runif(n())),
           "size_2"=quantile(HMP2.params.filt$size, probs=runif(n()))) %>%
    as.data.frame()
  
  # Create the cluster for parallel excecution
  cl <- makeCluster(6)
  registerDoSNOW(cl)
  
  # START INNER LOOP
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  df <- foreach(i=1:nrow(params_set), .combine="rbind",
                .packages=c("ToyModel"),
                .options.snow=list(progress=progress)) %dopar% {
                  
                  couple <-
                    ToyModel::toy_model(n=10^4, cor=params_set[i,"cor"], M=1,
                                        qdist=VGAM::qzinegbin,
                                        param=data.frame(
                                          "munb"=c(params_set[i,"munb_1"],params_set[i,"munb_2"]),
                                          "size"=c(params_set[i,"size_1"],params_set[i,"size_2"]),
                                          "pstr0"=c(params_set[i,"pstr0_1"],params_set[i,"pstr0_2"])
                                        ))
                  
                  random_HMP2_NorTA_i <- random_HMP2$NorTA
                  random_HMP2_NorTA_i[,25] <- couple$NorTA[,1]
                  random_HMP2_NorTA_i[,125] <- couple$NorTA[,2]
                  
                  cor_PCLR <- random_HMP2_NorTA_i %>% ToyModel::clr() %>% cor
                  
                  data.frame("iteration"=iter,
                             "cor_input"=params_set[i,"cor"],
                             "cor_normal"=couple$cor_normal[1,2],
                             "munb_1"=params_set[i,"munb_1"],
                             "munb_2"=params_set[i,"munb_2"],
                             "size_1"=params_set[i,"size_1"],
                             "size_2"=params_set[i,"size_2"],
                             "pstr0_1"=params_set[i,"pstr0_1"],
                             "pstr0_2"=params_set[i,"pstr0_2"],
                             "cor_NorTA_PCLR"=cor_PCLR[25,125])
                }
  # END INNER LOOP
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  stopCluster(cl)
  result <- bind_rows(result, df)
  # END OUTER LOOP
  #----------------------------------------------------------------------------#
}
saveRDS(result, "Sparsity_Effects_zinbin_rare_ecdf.rds")
