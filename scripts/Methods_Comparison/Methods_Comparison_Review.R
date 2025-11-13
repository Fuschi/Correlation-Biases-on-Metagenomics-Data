library(tidyverse)
library(SpiecEasi)
library(propr)
library(ggpubr)
library(grid)
library(ggplotify)
library(gridExtra)

# Load custom function
source("CLR.R")
source("TRIU.R")
source("LAYOUT_SIGNED.R")

# READ AND FILTER DATA
#------------------------------------------------------------------------------#

# Read HMP2
otu <- readRDS("../../data/otu_HMP2.rds")
meta <- readRDS("../../data/meta_HMP2.rds")
taxa <- readRDS("../../data/taxonomy.rds")

# Select samples belonging to 69-001 subject in health status
otu.69001.H <- otu[meta$SubjectID=="69-001" & meta$CL4_2=="Healthy", ]

# Remove rarest OTUs using prevalence and median of non-zero values
otu.filt <- otu.69001.H[, colSums(otu.69001.H>0)/nrow(otu.69001.H) >= .33]
otu.filt <- otu.filt[, apply(otu.filt,2,function(x) median(x[x>0])>=5)]
taxa.filt <- taxa[colnames(otu.filt), ]

# L1
res.L1  <- cor(otu.filt / rowSums(otu.filt), method="pearson")

# SparCC
res.cc  <- sparcc(otu.filt)$Cor
colnames(res.cc) <- rownames(res.cc) <- colnames(otu.filt)

# Rho
res.rho <- propr::propr(counts=otu.filt, metric="rho")@matrix

# Pearson+CLR
res.clr <- cor(CLR(otu.filt), method="pearson")

# PLOTS

# L1
p0 <- ggpubr::ggscatter(data.frame("PearsonCLR"=TRIU(res.clr),
                                   "PearsonL1"=TRIU(res.L1)),
                        x="PearsonCLR", y="PearsonL1",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), 
                   label.x=.45, label.y=-.25, size=6) +
  theme_bw() +
  xlab("Pearson+CLR") + ylab("Pearson+L1") +
  theme(plot.title = element_text(hjust = 0.5))


# SparCC
p1 <- ggpubr::ggscatter(data.frame("PearsonCLR"=TRIU(res.clr),
                                   "SparCC"=TRIU(res.cc)),
                        x="PearsonCLR", y="SparCC",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), 
                   label.x=.45, label.y=-.25, size=6) +
  theme_bw() +
  xlab("Pearson+CLR") +
  theme(plot.title = element_text(hjust = 0.5))

#Rho
p2 <- ggpubr::ggscatter(data.frame("PearsonCLR"=TRIU(res.clr),
                                   "Rho"=TRIU(res.rho)),
                        x="PearsonCLR", y="Rho",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label = after_stat(r.label)), label.x=.45, label.y=-.25, size=6) +
  theme_bw() +
  xlab("Pearson+CLR") +
  theme(plot.title = element_text(hjust = 0.5))

# Aggregate to phylum
#------------------------------------------------------------------------------#

phy <- otu.69001.H %>%
  as_tibble(rownames = "sample_id") %>%
  pivot_longer(-sample_id, names_to = "OTU", values_to = "abundance") %>%
  left_join(as_tibble(taxa.filt, rownames = "OTU"), by="OTU") %>%
  filter(!is.na(phylum)) %>%
  group_by(sample_id, phylum) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = phylum, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# L1
res.L1.phy  <- cor(phy / rowSums(phy), method="pearson")

# SparCC
res.cc.phy  <- sparcc(phy)$Cor
colnames(res.cc.phy) <- rownames(res.cc.phy) <- colnames(phy)

# Rho
res.rho.phy <- propr::propr(counts=phy, metric="rho")@matrix

# Pearson+CLR
res.clr.phy <- cor(CLR(phy), method="pearson")

# PLOTS

# L1
p0.phy <- ggpubr::ggscatter(data.frame("PearsonCLR"=TRIU(res.clr.phy),
                                       "PearsonL1"=TRIU(res.L1.phy)),
                        x="PearsonCLR", y="PearsonL1",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), 
                   label.x=-.3, label.y=-.75, size=6) +
  theme_bw() +
  xlab("Pearson+CLR") + ylab("Pearson+L1") 


# SparCC
p1.phy <- ggpubr::ggscatter(data.frame("PearsonCLR"=TRIU(res.clr.phy),
                                       "SparCC"=TRIU(res.cc.phy)),
                        x="PearsonCLR", y="SparCC",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), 
                   label.x=-0.7, label.y=.9, size=6) +
  theme_bw() +
  xlab("Pearson+CLR") 
#Rho
p2.phy <- ggpubr::ggscatter(data.frame("PearsonCLR"=TRIU(res.clr.phy),
                                   "Rho"=TRIU(res.rho.phy)),
                        x="PearsonCLR", y="Rho",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label = after_stat(r.label)), label.x=-.7, label.y=.15, size=6) +
  theme_bw() +
  xlab("Pearson+CLR") +
  theme(plot.title = element_text(hjust = 0.5))


# Combine plots

# Plot vuoto con etichetta "OTU level"
label_otu <- ggplot() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "OTU level", angle = 90, size = 6)

# Plot vuoto con etichetta "Phylum level"
label_phylum <- ggplot() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "Phylum level", angle = 90, size = 6)

combined_plots <- ggpubr::ggarrange(
  label_otu, p0, p2, p1,
  label_phylum, p0.phy, p2.phy, p1.phy,
  ncol = 4, nrow = 2,
  widths = c(0.15, 1, 1, 1), legend = "none",
  labels = c("", "A", "B", "C", "", "D", "E", "F")
)

png(filename = "../Plots/Methods_comparison_review.png", width = 13.8*600, height = 9*600, res = 600)
combined_plots
dev.off()
