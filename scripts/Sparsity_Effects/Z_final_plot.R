library(tidyverse)
library(igraph)
library(ggpubr)
library(magick)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#p_netw <- image_read("example_netw.png") %>% image_ggplot()
p_err <- readRDS("error.rds") + theme(legend.position="top")
p_couple <- readRDS("example_couple.rds")

#col1 <- ggarrange(p_netw, p_couple, ncol=1, labels=c("A","C"), heights=c(.4,.6))
#pall <- ggarrange(col1, p_err, ncol=2, labels=c("","B"))
pall <- ggarrange(p_err, p_couple, common.legend=T, legend="top")

png(filename="../Plots/sparsity.png", width=4800, height=2400, res=400)
pall
dev.off()
