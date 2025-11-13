library(tidyverse)
library(ggpubr)
library(cowplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df <- readRDS("Sparsity_Effects_htrlnorm_randzero.rds") %>%
  mutate("ERR_CLR"=abs(cor_normal-cor_NorTA_PCLR)) %>%
  mutate("phi_mean"=.5*(phi_1+phi_2)) %>%
  mutate("phi_max"=pmax(phi_1,phi_2)) %>%
  mutate("phi_min"=pmin(phi_1,phi_2))


if(any(df$ERR_CLR>1)) stop("Find Error greater than 1")

myPalette <- 
  RColorBrewer::brewer.pal(n=11, "Spectral") %>% rev() %>% 
  grDevices::colorRampPalette()

# CLR Effects
p.htrlnorm.min <- df %>%
  reframe(ERR_CLR_MEAN=mean(ERR_CLR), .by=c(cor_input,phi_min)) %>%
  ggplot(aes(x=cor_input, 
             y=phi_min,
             fill=ERR_CLR_MEAN)) +
  geom_tile() + theme_bw() + 
  scale_fill_gradientn(
    name="Absolute Error", 
    colours=myPalette(12), 
    values=c(seq(0,.5,by=.05),1),
    limits=c(0,1)) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(.1,.9,.1), expand=c(0,0))  +
  scale_x_continuous(breaks=seq(-.8,.8,.2), expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Correlation") + ylab("Zero %") 

p.htrlnorm.mean <- df %>%
  reframe(ERR_CLR_MEAN=mean(ERR_CLR), .by=c(cor_input,phi_mean)) %>% 
  filter((100*phi_mean)%%5==0) %>%
  ggplot(aes(x=cor_input, 
             y=factor(phi_mean),
             fill=ERR_CLR_MEAN)) +
  geom_tile() + theme_bw() + 
  scale_fill_gradientn(
    name="Absolute Error", 
    colours=myPalette(12), 
    values=c(seq(0,.5,by=.05),1),
    limits=c(0,1)) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_discrete(breaks=seq(.1,.9,.1), expand=c(0,0))  +
  scale_x_continuous(breaks=seq(-.8,.8,.2), expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Correlation") + ylab("Zero %") 

p.htrlnorm.max <- df %>%
  reframe(ERR_CLR_MEAN=mean(ERR_CLR), .by=c(cor_input,phi_max)) %>%
  ggplot(aes(x=cor_input, 
             y=phi_max,
             fill=ERR_CLR_MEAN)) +
  geom_tile() + theme_bw() + 
  scale_fill_gradientn(
    name="Absolute Error", 
    colours=myPalette(12), 
    values=c(seq(0,.5,by=.05),1),
    limits=c(0,1)) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(.1,.9,.1), expand=c(0,0))  +
  scale_x_continuous(breaks=seq(-.8,.8,.2), expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Correlation") + ylab("Zero %") 

png("intermediate/effects_htrlnorm_zerorand.png",width=3600, height=1500, res=300)
ggarrange(p.htrlnorm.min, p.htrlnorm.mean, p.htrlnorm.max,
          ncol=3, labels=c("Phi Min","Phi Mean","Phi Max"),
          common.legend=T, legend="bottom")
dev.off()




# Some Examples with htrlnorm
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Read HMP2
otu <- readRDS("../../data/otu_HMP2.rds")
meta <- readRDS("../../data/meta_HMP2.rds")
otu.69001.H <- otu[meta$SubjectID=="69-001" & meta$CL4_2=="Healthy", ]
otu.filt <- otu.69001.H[, colSums(otu.69001.H>0)/nrow(otu.69001.H) >= .33]
otu.filt <- otu.filt[, apply(otu.filt,2,function(x) median(x[x>0])>=5)]
#otu.filt <- round(otu.filt/rowSums(otu.filt) * min(rowSums(otu.filt)))
rm(otu, otu.69001.H, meta)
# Quantiles
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

example <- function(n, cor, meanlog_1, meanlog_2, sdlog_1, sdlog_2, phi_1, phi_2){
  
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
  
  random_HMP2 <- ToyModel::toy_model(n=n, cor=diag(200), M=1,
                                     qdist=HurdleTruncatedLogNormal::qhtrlnorm,
                                     param=params_random_HMP2)
  
  params_set <- data.frame("phi"=c(phi_1,phi_2),
                           "meanlog"=c(meanlog_1,meanlog_2),
                           "sdlog"=c(sdlog_1,sdlog_2)) %>%
    mutate(b=predict_max(meanlog))
  
  couple <- 
    ToyModel::toy_model(n=n, cor=cor, M=1,
                        qdist=HurdleTruncatedLogNormal::qhtrlnorm,
                        param=params_set)
  

  
  random_HMP2_NorTA <- random_HMP2$NorTA
  random_HMP2_NorTA[,1] <- couple$NorTA[,1]
  random_HMP2_NorTA[,2] <- couple$NorTA[,2]
  
  rand_pseudo <- 
    matrix(runif(length(random_HMP2_NorTA),min=.065, max=.65), 
           nrow=nrow(random_HMP2_NorTA), ncol=ncol(random_HMP2_NorTA))
  random_HMP2_NorTA <- random_HMP2_NorTA + (random_HMP2_NorTA==0)*rand_pseudo
  
  NorTA_CLR <- random_HMP2_NorTA %>% clr
  NorTA_CLR <- NorTA_CLR[,1:2]
  
  return(list("Normal"=couple$normal,
              "NorTA"=couple$NorTA,
              "NorTA_CLR"=NorTA_CLR))
}
plot_example <- function(example){
  
  p_norm <- example$Normal %>% as_tibble() %>%
    ggscatter(x="V1",y="V2", add="reg.line",
              add.params=list(color="blue", fill="lightgray"),
              conf.int = TRUE) +
    stat_cor(aes(label = ..r.label..), color="red") +
    ggtitle("Normal") + theme(axis.title=element_blank())
  
  p_NorTA <- example$NorTA %>% as_tibble() %>%
    ggscatter(x="V1",y="V2", add="reg.line",
              add.params=list(color="blue", fill="lightgray"),
              conf.int = TRUE) +
    stat_cor(aes(label = ..r.label..), color="red") +
    ggtitle("NorTA") + theme(axis.title=element_blank())
  
  p_NorTA_CLR <- example$NorTA_CLR %>% as_tibble() %>%
    ggscatter(x="V1",y="V2", add="reg.line",
              add.params=list(color="blue", fill="lightgray"),
              conf.int = TRUE) +
    stat_cor(aes(label = ..r.label..), color="red") +
    ggtitle("NorTA_CLR") + theme(axis.title=element_blank())
  
  return(list("norm"=p_norm, "NorTA"=p_NorTA, "NorTA_CLR"=p_NorTA_CLR))
}

e1 <- example(n=1000, cor=-.9, 
              meanlog_1=1, meanlog_2=3, sdlog_1=1, sdlog_2=1.3, 
              phi_1=.85, phi_2=.9) %>% plot_example()
pall1 <- ggarrange(e1$norm, e1$NorTA, e1$NorTA_CLR, ncol=3) 
pall1 <- annotate_figure(pall1,"n=1000, cor=-0.9, meanlog=[1,3], sdlog=[1,1.3], phi=[0.85,0.9]")

e2 <- example(n=1000, cor=.9, 
              meanlog_1=4, meanlog_2=2, sdlog_1=1.5, sdlog_2=.9, 
              phi_1=.9, phi_2=.75) %>% plot_example()
pall2 <- ggarrange(e2$norm, e2$NorTA, e2$NorTA_CLR, ncol=3) 
pall2 <- annotate_figure(pall2,"n=1000, cor=-0.9, meanlog=[4,2], sdlog=[1.5,0.9], phi=[0.9,0.75]")

e3 <- example(n=1000, cor=-.9, 
              meanlog_1=1.5, meanlog_2=5, sdlog_1=1, sdlog_2=1.5, 
              phi_1=.1, phi_2=.05) %>% plot_example()
pall3 <- ggarrange(e3$norm, e3$NorTA, e3$NorTA_CLR, ncol=3) 
pall3 <- annotate_figure(pall3,"n=1000, cor=-0.9, meanlog=[1.5,5], sdlog=[1,1.5], phi=[0.1,0.05]")

e4 <- example(n=1000, cor=.9, 
              meanlog_1=2.5, meanlog_2=5.5, sdlog_1=1.7, sdlog_2=2, 
              phi_1=.1, phi_2=.15) %>% plot_example()
pall4 <- ggarrange(e4$norm, e4$NorTA, e4$NorTA_CLR, ncol=3) 
pall4 <- annotate_figure(pall4,"n=1000, cor=0.9, meanlog=[2.5,5.5], sdlog=[1.7,2], phi=[0.1,0.15]")

pall <- plot_grid(pall1,pall2,pall3,pall4,nrow=4)

png("intermediate/examples_htrlnorm_zerorand.png", width=3600, height=4800, res=300)
pall
dev.off()