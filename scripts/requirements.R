# List of all the packages necessary to reproduce the analysis with the associated
# command to install them.

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
install.packages("tidyverse")

# gridExtra: arrange multiple grid-based plots on a page 
install.packages("gridExtra")

# SpiecEasi: package with spiec.easi and sparCC methods 
install_github("zdk123/SpiecEasi")

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

# lightweight package with necessary function to generate metagenomics data 
# based on the hurdle log-normal distribution and Normal To Anything method
install_github("https://github.com/Fuschi/ToyModel")
 