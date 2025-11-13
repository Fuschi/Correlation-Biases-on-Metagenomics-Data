library(tidyverse)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df <- readRDS("Compositional_Effects_2.rds")

df_sort <- tibble()
for(di in seq(5,200,by=5)){
  for(ei in seq(0.025,.975,by=.025)){
    
    sub_df <- df %>% filter(d==di) 
    rowi <- sub_df[which.min(abs(sub_df$pielou-ei)) ,]
    rowi <- c(rowi,"pielou_round"=ei)
    df_sort <- bind_rows(df_sort, rowi)
    
  }
}
rm(di,ei,sub_df,rowi)

df_sort <- df_sort %>%
  mutate(pielou_error=abs(pielou_round-pielou)) %>%
  mutate(pielou_error_logical=pielou_error<.005, .after=pielou_error) %>%
  mutate(LOG_ERR_L1=log10(ERR_L1)) %>%
  mutate_at("LOG_ERR_L1", \(x)(ifelse(x < -2, -2, x))) %>%
  mutate(LOG_ERR_CLR=log10(ERR_CLR)) %>%
  mutate_at("LOG_ERR_CLR", \(x)(ifelse(x < -2, -2, x)))

df_control <- df_sort %>%
  filter(pielou_error_logical==FALSE) 

myPalette <- RColorBrewer::brewer.pal(11, "Spectral") %>%
  rev() %>% grDevices::colorRampPalette()

df_sort %>% select(d,pielou_round,LOG_ERR_CLR) %>%
  pivot_wider(names_from=d, values_from=LOG_ERR_CLR) %>%
  column_to_rownames("pielou_round") %>% dim()


# L1 Effects
#png(filename="../outputs/L1_effects.png",width=3000,height=2400,res=600)
pL1 <- ggplot(df_sort, aes(x=d, y=pielou_round, fill=LOG_ERR_L1)) +
  geom_tile() + theme_bw() + 
  scale_fill_gradientn("MAE",limits=c(-2,0), colours=myPalette(11),
                       breaks=c(-0.1,-1,-2), labels=c(">1","0.1","< 0.01")) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(.05,.95,.05), expand=c(0,0))  +
  scale_x_continuous(breaks=seq(10,200,10), expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(expression(bar(P))) + xlab("D") +
  ggtitle("L1 Bias")
#pL1
#dev.off()

# CLR Effects
#png(filename="../outputs/CLR_effects.png",width=3000,height=2400,res=600)
pCLR <- ggplot(df_sort, aes(x=d, y=pielou_round, fill=LOG_ERR_CLR)) +
  geom_tile() + theme_bw() + 
  scale_fill_gradientn("MAE",limits=c(-2,0), colours=myPalette(11),
                       breaks=c(-0.1,-1,-2), labels=c(">1","0.1","< 0.01")) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_continuous(breaks=seq(.05,.95,.05), expand=c(0,0))  +
  scale_x_continuous(breaks=seq(10,200,10), expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(expression(bar(P))) + xlab("D") +
  ggtitle("CLR Bias")
#pCLR
#dev.off()


png(filename="../Plots/Normalization_Bias.png",width=6000, height=3000, res=600)
ggarrange(pL1, pCLR, labels=c("L1","CLR"), common.legend=T, 
          label.y=.125, label.x=c(0.03,-.02))
dev.off()


png(filename="../Plots/CLR_Compositional.png", width=1200, height=1200, res=300)
pCLR_Dim <- df_sort %>%
  reframe("ERR_CLR_meanD"=mean(ERR_CLR), .by=d) %>%
  ggplot(aes(x=d, y=ERR_CLR_meanD)) +
  geom_point(size=2.5) +
  geom_line(linewidth=1) +
  #coord_trans(y="log") +
  theme_bw() +
  #scale_y_log10() +
  scale_y_continuous(breaks=c(.01,.02,.05,.1,.15,.2)) +
  ylab("MAE") + xlab("D") +
  theme(plot.margin=unit(c(2,1,1,1),"cm"))
pCLR_Dim
dev.off()


## All plots merged
p_L1_CLR <- ggarrange(pL1, pCLR, common.legend=T, ncol=1)

p_all <- ggarrange(p_L1_CLR, 
                   pCLR_Dim + theme(axis.text=element_text(size=14),
                                    axis.title=element_text(size=16)), 
                   ncol=2, widths=c(.35,.65),
                   labels=c("A","B"), label.x=c(.05,.9))

png(filename="../Plots/Normalization_Bias_all.png", width=6000, height=4500, res=600)
p_all
dev.off()



df_percentiles <- df_sort %>%
  group_by(d) %>%
  summarise(
    ERR_CLR_mean = mean(ERR_CLR),
    ERR_CLR_p10 = quantile(ERR_CLR, 0.10),
    ERR_CLR_p90 = quantile(ERR_CLR, 0.90)
  )

pCLR_Dim <- ggplot(df_percentiles, aes(x=d, y=ERR_CLR_mean)) +
  geom_errorbar(aes(ymin = ERR_CLR_p10, ymax = ERR_CLR_p90), width = 2, color = "red") +
  #geom_line(linewidth=1.2, color="black") +
  geom_point(size=.5, color="black") +
  theme_bw() +
  scale_y_continuous(breaks=c(.01,.02,.05,.1,.15,.2)) +
  ylab("MAE") + xlab("D") +
  theme(plot.margin=unit(c(2,1,1,1),"cm"))

png(filename="../Plots/CLR_Compositional_percentiles.png", width=1200, height=1200, res=300)
pCLR_Dim
dev.off()