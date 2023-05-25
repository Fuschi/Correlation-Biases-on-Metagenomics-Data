library(tidyverse)

# Read HMP2 
otu <- readRDS("data/otu_HMP2.rds")

# Histogram of Depths
p1 <- ggplot(data=data.frame("depth"=log10(rowSums(otu))), aes(x=depth)) +
  geom_histogram(color="darkblue", fill="lightblue", bins=30) +
  theme_bw() + xlim(3,6) +
  xlab("log10(Depths)") +
  ggtitle("Samples Depth") +
  theme(plot.title = element_text(hjust=0.5))

# Presence heatmap
p2 <- ggplot(reshape2::melt(otu>0),aes(x=Var1,y=Var2,fill=value)) + 
  geom_tile() + theme_bw() + coord_flip() +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  labs(x="Taxa",y="Samples",title="Detection Image from Bacteria Matrix") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_discrete(name="Detection",labels=c("Absence  (=0)", "Presence (>0)")) +
  theme(legend.key.size = unit(.4, 'cm')) +
  theme(legend.position = c(.835,.89)) +
  ggtitle("Counts Sparsity") +
  theme(plot.title = element_text(hjust=0.5))

png(filename="outputs/NGS_Features.png",width=2400,height=1200, res=300)
gridExtra::grid.arrange(p1,p2,ncol=2)
dev.off()