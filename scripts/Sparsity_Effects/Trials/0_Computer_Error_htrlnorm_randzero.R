library(tidyverse)
library(doSNOW)
library(progress)
library(foreach)
library(ToyModel)

setwd("~/Scrivania/Correlation-Biases-on-Metagenomics-Data/scripts/Sparsity_Effects/")

# Read HMP2
otu <- readRDS("../../data/otu_HMP2.rds")
meta <- readRDS("../../data/meta_HMP2.rds")

# Select samples belonging to 69-001 subject in health status
otu.69001.H <- otu[meta$SubjectID=="69-001" & meta$CL4_2=="Healthy", ]

# Remove rarest OTUs using prevalence and median of non-zero values
otu.filt <- otu.69001.H[, colSums(otu.69001.H>0)/nrow(otu.69001.H) >= .33]
otu.filt <- otu.filt[, apply(otu.filt,2,function(x) median(x[x>0])>=5)]
#otu.filt <- round(otu.filt/rowSums(otu.filt) * min(rowSums(otu.filt)))
rm(otu, otu.69001.H, meta)

# Fit ZINB params from real data for each filtered OTUs and
# save the first and ninen-th decil of the distribution
HMP2.quantile.params <- apply(otu.filt,2,function(x){
  HurdleTruncatedLogNormal::mle.htrlnorm(x)$estimate
}) %>% apply(1,function(x) quantile(x,probs=c(.1,.9)))

# Fitting a linear model between the mean and max in log scale for each OTU
df_mean_max <- tibble(
  "y"=apply(log(otu.filt+1),2,max),
  "x"=apply(log(otu.filt+1),2,mean)
)
model <- lm(y~x, data=df_mean_max)

# Function to predict new maximum values from the means
predict_max <- function(new_xs){
  
  # Predict Mean Values and Standard Errors
  new_data <- data.frame(x=new_xs)
  predicted_values <- predict(model, new_data, interval = "none")
  
  # Calculate standard error of prediction
  residuals_variance <- sum(residuals(model)^2) / model$df.residual
  n <- length(model$model$x)  # Number of observations in original data
  x_bar <- mean(model$model$x)  # Mean of original independent variable
  
  # Leverage for each new_x (hii = 1/n + (xi - x̄)^2 / Σ(xi - x̄)^2)
  leverages <- 1/n + ((new_data$x - x_bar)^2 / sum((model$model$x - x_bar)^2))
  
  # Standard error of prediction for each new_x
  std_error_prediction <- sqrt(residuals_variance * (1 + leverages))
  
  
  # Generate random values from normal distributions with these means and variances
  simulated_ys <- rnorm(nrow(new_data), mean = predicted_values, sd = std_error_prediction)
  return(simulated_ys)
}

# Create the progression bar
nIteration <- 100
pb <- progress_bar$new(
  format = "[:bar] :elapsed | eta: :eta",
  total = nIteration*7600,
  width = 60)
progress <- function(n){pb$tick()}

set.seed(42)
result <- data.frame()
for(iter in 1:nIteration){
  # START OUTER LOOP
  #----------------------------------------------------------------------------#
  
  # GENERATION OF THE UNDERLINE DATASET
  params_random_HMP2 <- data.frame(
    "meanlog"=runif(n=200,
                    min=HMP2.quantile.params["10%","meanlog"],
                    max=HMP2.quantile.params["90%","meanlog"]),
    "sdlog"=runif(n=200,
                  min=HMP2.quantile.params["10%","sdlog"],
                  max=HMP2.quantile.params["90%","sdlog"]),
    "phi"=runif(n=200,
                min=HMP2.quantile.params["10%","phi"],
                max=HMP2.quantile.params["90%","phi"]))
  
  params_random_HMP2$b <- predict_max(params_random_HMP2$meanlog)
  

  random_HMP2 <- ToyModel::toy_model(n=10^4, cor=diag(200), M=1,
                                     qdist=HurdleTruncatedLogNormal::qhtrlnorm,
                                     param=params_random_HMP2)
  random_cor0_HMP2 <- random_HMP2$cor_normal
  

  # CREATE PARAMETERS OF THE SUPERVISED COUPLE OF VARIABLES
  params_set <- expand_grid(
    "phi_1"=seq(0,.95,by=.05),
    "phi_2"=seq(0,.95,by=.05),
    "cor"=seq(-.9,.9,by=.1)
  ) %>%
    mutate("meanlog_1"=runif(n=n(),
                          min=HMP2.quantile.params["10%","meanlog"],
                          max=HMP2.quantile.params["90%","meanlog"]),
           "meanlog_2"=runif(n=n(),
                          min=HMP2.quantile.params["10%","meanlog"],
                          max=HMP2.quantile.params["90%","meanlog"]),
           "sdlog_1"=runif(n=n(),
                          min=HMP2.quantile.params["10%","sdlog"],
                          max=HMP2.quantile.params["90%","sdlog"]),
           "sdlog_2"=runif(n=n(),
                          min=HMP2.quantile.params["10%","sdlog"],
                          max=HMP2.quantile.params["90%","sdlog"])) %>%
    as.data.frame() %>%
    mutate(b_1=predict_max(meanlog_1),
           b_2=predict_max(meanlog_2))
  
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
                                        qdist=HurdleTruncatedLogNormal::qhtrlnorm,
                                        param=data.frame(
                                          "meanlog"=c(params_set[i,"meanlog_1"],params_set[i,"meanlog_2"]),
                                          "sdlog"=c(params_set[i,"sdlog_1"],params_set[i,"sdlog_2"]),
                                          "phi"=c(params_set[i,"phi_1"],params_set[i,"phi_2"])
                                        ))
                  
                  random_HMP2_NorTA_i <- random_HMP2$NorTA
                  random_HMP2_NorTA_i[,25] <- couple$NorTA[,1]
                  random_HMP2_NorTA_i[,125] <- couple$NorTA[,2]
                  
                  # Replace Zero with random pseudo-values between 0.1*dl and dl
                  rand_pseudo <- 
                    matrix(runif(length(random_HMP2_NorTA_i),min=.065, max=.65), 
                           nrow=nrow(otu.filt), ncol=ncol(otu.filt))
                  random_HMP2_NorTA_i <- random_HMP2_NorTA_i + (random_HMP2_NorTA_i==0)*rand_pseudo
                  
                  cor_PCLR <- random_HMP2_NorTA_i %>% ToyModel::clr() %>% cor
                  
                  data.frame("iteration"=iter,
                             "cor_input"=params_set[i,"cor"],
                             "cor_normal"=couple$cor_normal[1,2],
                             "meanlog_1"=params_set[i,"meanlog_1"],
                             "meanlog_2"=params_set[i,"meanlog_2"],
                             "sdlog_1"=params_set[i,"sdlog_1"],
                             "sdlog_2"=params_set[i,"sdlog_2"],
                             "b_1"=params_set[i,"b_1"],
                             "b_2"=params_set[i,"b_2"],
                             "phi_1"=params_set[i,"phi_1"],
                             "phi_2"=params_set[i,"phi_2"],
                             "cor_NorTA_PCLR"=cor_PCLR[25,125])
                }
  # END INNER LOOP
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  stopCluster(cl)
  result <- bind_rows(result, df)
  # END OUTER LOOP
  #----------------------------------------------------------------------------#
}
saveRDS(result, "Sparsity_Effects_htrlnorm_rand_pseudo.rds")
