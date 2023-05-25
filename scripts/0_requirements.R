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
