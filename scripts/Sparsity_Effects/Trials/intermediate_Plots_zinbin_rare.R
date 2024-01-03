library(tidyverse)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df <- readRDS("Sparsity_Effects_zinbin_rare.rds") %>%
  mutate("ERR_CLR"=abs(cor_normal-cor_NorTA_PCLR)) %>%
  mutate("pstr0_mean"=.5*(pstr0_1+pstr0_2)) %>%
  mutate("pstr0_max"=pmax(pstr0_1,pstr0_2)) %>%
  mutate("pstr0_min"=pmin(pstr0_1,pstr0_2))

if(any(df$ERR_CLR>1)) stop("Find Error greater than 1")

myPalette <- 
  RColorBrewer::brewer.pal(n=11, "Spectral") %>% rev() %>% 
  grDevices::colorRampPalette()

# CLR Effects
p.zinb.rare.min <- df %>%
  reframe(ERR_CLR_MEAN=mean(ERR_CLR), .by=c(cor_input,pstr0_min)) %>%
  ggplot(aes(x=cor_input, 
             y=pstr0_min,
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

p.zinb.rare.mean <- df %>%
  reframe(ERR_CLR_MEAN=mean(ERR_CLR), .by=c(cor_input,pstr0_mean)) %>% 
  filter((100*pstr0_mean)%%5==0) %>%
  ggplot(aes(x=cor_input, 
             y=factor(pstr0_mean),
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

p.zinb.rare.max <- df %>%
  reframe(ERR_CLR_MEAN=mean(ERR_CLR), .by=c(cor_input,pstr0_max)) %>%
  ggplot(aes(x=cor_input, 
             y=pstr0_max,
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

png("intermediate/effects_zinb_rare.png",width=3600, height=1500, res=300)
ggarrange(p.zinb.rare.min, p.zinb.rare.mean, p.zinb.rare.max,
          ncol=3, labels=c("pstr0 Min","pstr0 Mean","pstr0 Max"),
          common.legend=T, legend="bottom")
dev.off()

# saveRDS(p.zinb.rare.min, "intermediate/p_zinb_rare_min.rds")
# saveRDS(p.zinb.rare.mean, "intermediate/p_zinb_rare_mean.rds")
# saveRDS(p.zinb.rare.max, "intermediate/p_zinb_rare_max.rds")

# Some Examples with zinb
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Read HMP2
otu <- readRDS("../../data/otu_HMP2.rds")
meta <- readRDS("../../data/meta_HMP2.rds")
otu.69001.H <- otu[meta$SubjectID=="69-001" & meta$CL4_2=="Healthy", ]
otu.filt <- otu.69001.H[, colSums(otu.69001.H>0)/nrow(otu.69001.H) >= .33]
otu.filt <- otu.filt[, apply(otu.filt,2,function(x) median(x[x>0])>=5)]
otu.filt <- round(otu.filt/rowSums(otu.filt) * min(rowSums(otu.filt)))
rm(otu, otu.69001.H, meta)
# Quantiles
HMP2.params <- apply(otu.filt,2,function(x){
  SpiecEasi::fitdistr(as.numeric(x),"zinegbin")$par
}) 

HMP2.quantile.params <- HMP2.params %>%
  apply(1,function(x) quantile(x,probs=c(.1,.9)))

example <- function(n, cor, munb_1, munb_2, size_1, size_2, pstr0_1, pstr0_2){
  
  params_random_HMP2 <- data.frame(
    "munb"=runif(n=200,
                    min=HMP2.quantile.params["10%","munb"],
                    max=HMP2.quantile.params["90%","munb"]),
    "size"=runif(n=200,
                  min=HMP2.quantile.params["10%","size"],
                  max=HMP2.quantile.params["90%","size"]),
    "pstr0"=runif(n=200,
                min=HMP2.quantile.params["10%","pstr0"],
                max=HMP2.quantile.params["90%","pstr0"]))
  
  random_HMP2 <- ToyModel::toy_model(n=n, cor=diag(200), M=1,
                                     qdist=VGAM::qzinegbin,
                                     param=params_random_HMP2)
  
  params_set <- data.frame("pstr0"=c(pstr0_1,pstr0_2),
                           "munb"=c(munb_1,munb_2),
                           "size"=c(size_1,size_2)) 
  
  couple <- 
    ToyModel::toy_model(n=n, cor=cor, M=1,
                        qdist=VGAM::qzinegbin,
                        param=params_set)
  
  
  
  random_HMP2_NorTA <- random_HMP2$NorTA
  random_HMP2_NorTA[,1] <- couple$NorTA[,1]
  random_HMP2_NorTA[,2] <- couple$NorTA[,2]
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
              munb_1=20, munb_2=50, size_1=25, size_2=3, 
              pstr0_1=.85, pstr0_2=.9) %>% plot_example()
pall1 <- ggarrange(e1$norm, e1$NorTA, e1$NorTA_CLR, ncol=3) 
pall1 <- annotate_figure(pall1,"n=1000, cor=-0.9, munb=[20,50], size=[25,3], pstr0=[0.85,0.9]")

e2 <- example(n=1000, cor=.9, 
              munb_1=23, munb_2=15, size_1=3.5, size_2=5, 
              pstr0_1=.9, pstr0_2=.75) %>% plot_example()
pall2 <- ggarrange(e2$norm, e2$NorTA, e2$NorTA_CLR, ncol=3) 
pall2 <- annotate_figure(pall2,"n=1000, cor=0.9, munb=[23,15], size=[3.5,5], pstr0=[0.9,0.75]")

e3 <- example(n=1000, cor=-.9, 
              munb_1=20, munb_2=45, size_1=5.5, size_2=2.5, 
              pstr0_1=.1, pstr0_2=.05) %>% plot_example()
pall3 <- ggarrange(e3$norm, e3$NorTA, e3$NorTA_CLR, ncol=3) 
pall3 <- annotate_figure(pall3,"n=1000, cor=-0.9, munb=[20,45], size=[5.5,2.5], pstr0=[0.1,0.05]")

e4 <- example(n=1000, cor=.9, 
              munb_1=15, munb_2=11, size_1=1.7, size_2=0.7, 
              pstr0_1=.1, pstr0_2=.15) %>% plot_example()
pall4 <- ggarrange(e4$norm, e4$NorTA, e4$NorTA_CLR, ncol=3) 
pall4 <- annotate_figure(pall4,"n=1000, cor=0.9, munb=[15,11], size=[1.7,0.7], pstr0=[0.1,0.15]")

e5 <- example(n=1000, cor=.9, 
              munb_1=24, munb_2=15, size_1=3.5, size_2=5, 
              pstr0_1=0, pstr0_2=0) %>% plot_example()
pall5 <- ggarrange(e5$norm, e5$NorTA, e5$NorTA_CLR, ncol=3) 
pall5 <- annotate_figure(pall5,"n=1000, cor=0.9, munb=[24,15], size=[3.5,5], pstr0=[0,0]")

e6 <- example(n=1000, cor=-.9, 
              munb_1=73, munb_2=24, size_1=0.5, size_2=1.9, 
              pstr0_1=0, pstr0_2=0) %>% plot_example()
pall6 <- ggarrange(e6$norm, e6$NorTA, e6$NorTA_CLR, ncol=3) 
pall6 <- annotate_figure(pall6,"n=1000, cor=0.9, munb=[73,24], size=[0.5,1.9], pstr0=[0,0]")

pall <- plot_grid(pall1,pall2,pall3,pall4,pall5,pall6,nrow=6)

png("intermediate/examples_zinb_rare.png", width=3600, height=7200, res=300)
pall
dev.off()
