library(tidyverse)
library(ggpubr)
library(zCompositions)

# Load custom function
source("../Methods_Comparison/CLR.R")

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

# Different Zero Strategies
#------------------------------------------------------------------------------#
cor(CLR(as.matrix(zCompositions::cmultRepl(otu.filt, method="CZM", label=0)))) -> cor.czm
cor(CLR(as.matrix(zCompositions::cmultRepl(otu.filt, method="GBM", label=0)))) -> cor.gbm
cor(CLR(as.matrix(zCompositions::cmultRepl(otu.filt, method="BL", label=0)))) -> cor.bl
otu.filt.65 <- otu.filt; otu.filt.65[otu.filt.65==0] <- .65; cor.65 <- cor(CLR(otu.filt.65))


# Nomi delle variabili (es. OTUs)
otu_names <- colnames(otu.filt)
idx <- which(upper.tri(cor.czm), arr.ind = TRUE)

# Crea un tibble tidy con i valori
cor_df <- tibble(
  var1 = otu_names[idx[, 1]],
  var2 = otu_names[idx[, 2]],
  CZM = cor.czm[upper.tri(cor.czm)],
  GBM = cor.gbm[upper.tri(cor.gbm)],
  BL  = cor.bl[upper.tri(cor.bl)],
  PC65 = cor.65[upper.tri(cor.65)]
)

cor_long <- cor_df |>
  pivot_longer(cols = c(CZM, GBM, BL, PC65), names_to = "Method", values_to = "Correlation")


tbl <- cor_long |>
  group_by(var1, var2) |>
  summarise(max_abs_diff = max(Correlation) - min(Correlation),
            mean_corr = mean(Correlation), .groups = "drop") %>%
  mutate(higher = ifelse(max_abs_diff > .1, TRUE, FALSE))


png(width = 1600, height = 1200, res = 300, filename = "../Plots/correlation_differences_between_zeroRepl.png")
tbl |>
  ggplot(aes(x = max_abs_diff)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "darkblue") +
  theme_bw() +
  xlab("Maximum Absolute Difference in Correlation\nbetween Zero-Replacement Strategies\n(CZM, GBM, BL, 65% detection threshold)")
dev.off()



