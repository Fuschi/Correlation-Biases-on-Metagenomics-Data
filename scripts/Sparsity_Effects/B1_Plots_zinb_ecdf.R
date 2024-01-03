library(tidyverse)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df <- readRDS("Sparsity_Effects_zinbin_rare_ecdf_2.rds") %>%
  dplyr::select(-pstr0_2) %>% rename(pstr0=pstr0_1) %>%
  mutate("ERR_CLR"=abs(cor_normal-cor_NorTA_PCLR))

if(any(df$ERR_CLR>1)) stop("Find Error greater than 1")

myPalette <- 
  c(RColorBrewer::brewer.pal(n=11, "Spectral")) %>% rev() %>% c(.,"#000000") %>%
  grDevices::colorRampPalette()

p <- df %>%
  reframe(MEAN_ERR_CLR=mean(ERR_CLR), .by=c(cor_input, pstr0)) %>%
  ggplot(aes(x=cor_input, y=pstr0, fill=MEAN_ERR_CLR)) +
  geom_tile() + theme_bw() +
  scale_fill_gradientn("MAE", colours=myPalette(12),
                       values=c(seq(0,.5,by=.05),1),
                       limits=c(0,1),
                       labels=c("0","0.25","0.5","0.75","1")) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(0,.9,.1), expand=c(0,0))  +
  scale_x_continuous(breaks=seq(-.8,.8,.2), expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Correlation") + ylab("Zero %") 

saveRDS(p, "error.rds")
