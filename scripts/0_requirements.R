# List of all the packages necessary to reproduce the analysis with the associated
# command to install them.

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
install.packages("tidyverse")

# reshape2: converts matrices to long format for sorting with the melt function
install.packages("reshape2")

# gridExtra: arrange multiple grid-based plots on a page 
install.packages("gridExtra")

# SpiecEasi: package with spiec.easi and sparCC methods 
install_github("zdk123/SpiecEasi")

# psych: computes p-value corrections on correlations measures
install.packages('psych')

# propr: package with proportionality rho method
install_github("https://github.com/tpq/propr")

# ggpubr: nice plots based on ggplot2
install.packages("ggpubr")

# ggVennDiagram: Venn diagram based on ggplot2
install.packages("ggVennDiagram")

# qualpalr: generate distinct qualitative color palette 
# (useful to associate vertex color to taxonomy)
install.packages("qualpalr")

# cowplot: draw ggplot2 in new figures
install.packages("cowplot")

# vegan: used to elaborates shannon entropy
install.packages("vegan")

# mvtnorm: Multivariate Normal utilities
install.packages("mvtnorm")
 