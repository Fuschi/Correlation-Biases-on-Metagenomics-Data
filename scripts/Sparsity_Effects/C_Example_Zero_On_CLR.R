library(tidyverse)
library(ggpubr)

setwd("~/Scrivania/Correlation-Biases-on-Metagenomics-Data/scripts/Sparsity_Effects/")

# Read HMP2
otu <- readRDS("../../data/otu_HMP2.rds")
meta <- readRDS("../../data/meta_HMP2.rds")

otu.filt <- otu[, colSums(otu>0)/nrow(otu)>=.25]
otu.filt <- otu.filt[, apply(otu.filt,2,\(x)median(x[x>0])>=5)]

otu.filt.CLR <- ToyModel::clr(otu.filt)
PCLR <- cor(otu.filt.CLR)

otu_prev <- setNames(colSums(otu.filt>0)/nrow(otu.filt), colnames(otu.filt))


info <- data.frame(PCLR) %>% rownames_to_column("OTU_I") %>%
  pivot_longer(!OTU_I, names_to="OTU_J", values_to="cor") %>%
  filter(OTU_I>OTU_J) %>%
  mutate(prev_I=otu_prev[OTU_I],
         prev_J=otu_prev[OTU_J]) %>%
  mutate(zero_I=100*round(1-prev_I,2),
         zero_J=100*round(1-prev_J,2))

info_filt <- info %>%
  filter(prev_I<=.5, prev_J<=.5, abs(cor)>=.4)

idx.max <- c("OTU_269","OTU_97")
idx.min <- c("OTU_269","OTU_313")

detection <- otu.filt %>% as_tibble() %>%
  select(all_of(idx.min)) %>%
  mutate(detection=case_when(
    OTU_269==0 & OTU_313==0 ~ "Both are 0",
    OTU_269==0 ~ "OTU_269 is 0",
    OTU_313==0 ~ "OTU_313 is 0",
    OTU_269>0 & OTU_313>0 ~ "Both are >0", # Adjusted condition for clarity
  )) %>%
  select(detection)

p.min <- otu.filt %>% as_tibble() %>%
  select(all_of(idx.min)) %>%
  cbind(detection) %>%
  ggscatter(x="OTU_269", y="OTU_313", add="reg.line", color="detection",
            size=1.5, palette=c("Both are 0"=rgb(.6,0,0,.4), 
                                "OTU_269 is 0"=rgb(0,.6,.3,.4), 
                                "OTU_313 is 0"=rgb(0,.3,.6,.4), 
                                "Both are >0"=rgb(0,0,0,.6)),
            add.params=list(color="#800080", fill="lightgray"),
            conf.int=TRUE) +
  stat_cor(aes(label = after_stat(r.label)), color="red", label.x.npc=0) +
  xlab("OTU 269 (Zero ~73%)") +
  ylab("OTU 313 (Zero ~55%)") +
  theme(axis.title=element_text(size=12))

p.min.clr <- otu.filt.CLR %>% as_tibble %>%
  dplyr::select(all_of(idx.min)) %>%
  cbind(detection) %>%
  ggscatter(x="OTU_269", y="OTU_313", add="reg.line", color="detection",
            size=1.5, palette=c("Both are 0"=rgb(.6,0,0,.4), 
                                "OTU_269 is 0"=rgb(0,.6,.3,.4), 
                                "OTU_313 is 0"=rgb(0,.3,.6,.4), 
                                "Both are >0"=rgb(0,0,0,.6)),
            add.params=list(color="#800080", fill="lightgray"),
            conf.int=TRUE) +
  stat_cor(aes(label = after_stat(r.label)), color="red", label.x.npc=0) +
  xlab("OTU 269 (Zero ~73%)") +
  ylab("OTU 313 (Zero ~55%)") +
  theme(axis.title=element_text(size=12))


detection <- otu.filt %>% as_tibble() %>%
  select(all_of(idx.max)) %>%
  mutate(detection=case_when(
    OTU_269==0 & OTU_97==0 ~ "Both are 0",
    OTU_269==0 ~ "OTU_269 is 0",
    OTU_97==0 ~ "OTU_97 is 0",
    OTU_269>0 & OTU_97>0 ~ "Both are >0", 
  )) %>%
  select(detection)

p.min <- otu.filt %>% as_tibble() %>%
  select(all_of(idx.min)) %>%
  cbind(detection) %>%
  ggscatter(x="OTU_269", y="OTU_313", add="reg.line", color="detection",
            size=1.5, palette=c("Both are 0"=rgb(.6,0,0,.4), 
                                "OTU_269 is 0"=rgb(0,.6,.3,.4), 
                                "OTU_313 is 0"=rgb(0,.3,.6,.4), 
                                "Both are >0"=rgb(0,0,0,.6)),
            add.params=list(color="#800080", fill="lightgray"),
            conf.int=TRUE) +
  stat_cor(aes(label = after_stat(r.label)), color="red", label.x.npc=0) +
  xlab("OTU 269 (Zero ~73%)") +
  ylab("OTU 313 (Zero ~55%)") +
  theme(axis.title=element_text(size=12))

p.min.clr <- otu.filt.CLR %>% as_tibble %>%
  dplyr::select(all_of(idx.min)) %>%
  cbind(detection) %>%
  ggscatter(x="OTU_269", y="OTU_313", add="reg.line", color="detection",
            size=1.5, palette=c("Both are 0"=rgb(.6,0,0,.4), 
                                "OTU_269 is 0"=rgb(0,.6,.3,.4), 
                                "OTU_313 is 0"=rgb(0,.3,.6,.4), 
                                "Both are >0"=rgb(0,0,0,.6)),
            add.params=list(color="#800080", fill="lightgray"),
            conf.int=TRUE) +
  stat_cor(aes(label = after_stat(r.label)), color="red", label.x.npc=0) +
  xlab("OTU 269 (Zero ~73%)") +
  ylab("OTU 313 (Zero ~55%)") +
  theme(axis.title=element_text(size=12))

p.max <- otu.filt %>% as_tibble %>%
  dplyr::select(all_of(idx.max)) %>%
  ggscatter(x="OTU_269",y="OTU_97", add="reg.line",
            size=1.5, color=rgb(0,0,0,.4),
            add.params=list(color="blue", fill="lightgray"),
            conf.int=TRUE) +
  stat_cor(aes(label = after_stat(r.label)), color="red", label.x.npc=.3) +
  xlab("OTU 269 (Zero ~73%)") +
  ylab("OTU 97 (Zero ~56%)") +
  theme(axis.title=element_text(size=6))

p.max.clr <- otu.filt.CLR %>% as_tibble %>%
  dplyr::select(all_of(idx.max)) %>%
  ggscatter(x="OTU_269",y="OTU_97", add="reg.line",
            size=1.5, color=rgb(0,0,0,.4),
            add.params=list(color="blue", fill="lightgray"),
            conf.int=TRUE) +
  stat_cor(aes(label = after_stat(r.label)), color="red", label.x.npc=.3) +
  xlab("OTU 269 (Zero ~73%)") +
  ylab("OTU 97 (Zero ~56%)") +
  theme(axis.title=element_text(size=6))


pall <- ggarrange(p.min, p.min.clr, p.max, p.max.clr, ncol=2, nrow=2)

saveRDS(pall,"example_couple.rds")
