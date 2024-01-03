library(igraph)
library(extrafont)

font_import()
loadfonts(device = "postscript")

setwd("~/Scrivania/Correlation-Biases-on-Metagenomics-Data/scripts/Sparsity_Effects/")

adj <- matrix(0, nrow=40, ncol=40)
adj[5,25] <- adj[25,5] <- 1

g <- graph_from_adjacency_matrix(adj, mode="undirected")
E(g)$color <- "darkgray"
E(g)$width <- 2 
E(g)$lty <- 2

set.seed(4)
V(g)$size <- runif(n=40, min=3, max=12.5)

png(filename="example_netw.png", width=1200, height=1200, res=300)
par(mar=c(0,0,0,0))
plot(g, layout=layout_in_circle(g), vertex.label=NA)
text(0.925,0.7,"I",family="Times New Roman",cex=2)
text(-0.925,-0.7,"J",family="Times New Roman",cex=2)
dev.off()



