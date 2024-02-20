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


# COMPUTES ALL METHODS AND PREPARE RESULTS FOR PLOTS
#------------------------------------------------------------------------------#

# Spiec-Easi glasso
res.gl <- spiec.easi(data=otu.filt, method='glasso', lambda.max=.75, lambda.min.ratio=.5,
                   pulsar.params=list(ncores=6, thresh=0.05))
adj.gl <- cov2cor(as.matrix(getOptCov(res.gl)))
adj.gl <- adj.gl * as.matrix(getRefit(res.gl))
colnames(res.gl$est$data) -> colnames(adj.gl) -> rownames(adj.gl)
adj.gl[adj.gl>0] <- 1
adj.gl[adj.gl<0] <- -1

# SparCC
res.cc  <- sparcc(otu.filt)$Cor
colnames(res.cc) <- rownames(res.cc) <- colnames(otu.filt)

# Rho
res.rho <- propr::propr(counts=otu.filt, metric="rho")@matrix

# Pearson+CLR
res.clr <- cor(CLR(otu.filt), method="pearson")
diag(res.clr) <- 0
p.adjust <- psych::corr.p(r=res.clr, n=nrow(otu.filt),
                  adjust="bonferroni", ci=FALSE)$p
p.adjust[lower.tri(p.adjust)] <- t(p.adjust)[lower.tri(p.adjust)]
diag(p.adjust) <- 1
adj.clr <- res.clr*(p.adjust<=.05)

# Re-arrange data for spiec-easi and Pearson+CLR histogram comparison
#------------------------------------------------------------------------------#
#res.mb.clr <- TRIU(res.clr*abs(adj.mb)) %>% .[.!=0]
res.gl.clr <- TRIU(res.clr*abs(adj.gl)) %>% .[.!=0]

# df.mb <- tibble("method"=rep("clr",length(TRIU(res.clr))),"value"=TRIU(res.clr)) %>%
#   rbind(tibble("method"=rep("mb",length(res.mb.clr)),"value"=res.mb.clr))

df.gl <- tibble("method"=rep("clr",length(TRIU(res.clr))),"value"=TRIU(res.clr)) %>%
  rbind(tibble("method"=rep("gl",length(res.gl.clr)),"value"=res.gl.clr))



# PLOTS
#------------------------------------------------------------------------------#

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


#GLASSO
p4.all<-ggplot(df.gl, aes(x=value, fill=method, color=method)) +
  geom_histogram(position="identity", alpha=.5, breaks=seq(-1,1,by=.1)) +
  geom_hline(yintercept=200, linetype="twodash", color="black", linewidth=1) +
  xlim(c(-1,1)) +
  theme_bw() +
  scale_color_manual(values=c("steelblue1", "forestgreen"))+
  scale_fill_manual(values=c("steelblue1", "forestgreen")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="right",
        legend.direction="vertical")

p4.zoom <- ggplot(df.gl, aes(x=value, fill=method, color=method)) +
  geom_histogram(position="identity", alpha=.5, breaks=seq(-1,1,by=.1)) +
  coord_cartesian(ylim=c(0, 200)) +
  theme_bw() +
  scale_color_manual(values=c("steelblue1", "forestgreen"))+
  scale_fill_manual(values=c("steelblue1", "forestgreen")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme_void() +
  theme(legend.position="none") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = +Inf), 
            col="black", alpha=0, linewidth=1)

p4 <- p4.all +
  annotation_custom(ggplotGrob(p4.zoom), xmin=-1.1, xmax=-.3, 
                    ymin=1400, ymax = 1900) 



# GRAPHS COMPARISON
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

g.gl <- igraph::graph_from_adjacency_matrix(adj.gl, mode="undirected",
                                            weighted=TRUE)

g.clr <- igraph::graph_from_adjacency_matrix(adj.clr, mode="undirected",
                                             weighted=TRUE)


# Venn Diagramm GGplot
p.venn <- ggVennDiagram::ggVennDiagram(x=list(
  "GLASSO"=paste(igraph::as_edgelist(g.gl)[,1],"-",igraph::as_edgelist(g.gl)[,2],sep=""),
  "CLR"=paste(igraph::as_edgelist(g.clr)[,1],"-",igraph::as_edgelist(g.clr)[,2],sep="")
), 
label_alpha=0) +
  ggplot2::scale_fill_gradient("Shared \n Links",low="white",high = "red") +
  ggplot2::scale_color_manual(values = c("gray20","gray20","gray20")) +
  theme(legend.position="right")

# Plot Networks
#------------------------------------------------------------------------------#

# Generate colormap for taxonomy at family level
set.seed(42)
colpal <- rownames(qualpalr::qualpal(n=length(unique(taxa.filt[,"family"])))$RGB)
names(colpal) <- unique(taxa.filt[,"family"])

# Vertex size proportional to clr abundances
vertex.size <- colMeans(CLR(otu.filt) - min(CLR(otu.filt)))

# CLR
set.seed(42)
p.graph.CLR <- as.grob(~plot(g.clr, vertex.label=NA,
                             vertex.color=colpal[taxa.filt[,"family"]],
                             vertex.size=vertex.size,
                             edge.color=ifelse(igraph::E(g.clr)$weight>0,rgb(0,0,1),rgb(1,0,0)),
                             edge.width=.5,
                             layout=LAYOUT_SIGNED(g.clr)))

# GLASSO
set.seed(42)
p.graph.glasso <- as.grob(~plot(g.gl, vertex.label=NA,
     vertex.color=colpal[taxa.filt[,"family"]],
     vertex.size=vertex.size,
     edge.color=ifelse(igraph::E(g.gl)$weight>0,rgb(0,0,1),rgb(1,0,0)),
     edge.width=.5,
     layout=LAYOUT_SIGNED(g.clr)))

f.leg <- function(){
  plot.new()
  par(mar=c(0,0,0,0))
  legend(x="center",legend=names(colpal), fill=colpal, ncol=1, cex=.6,
         xpd = TRUE, y.intersp=2, lwd=5)
}

p.leg.all <- grid.arrange(list(p.venn, f.leg()), ncol=1, heights=c(.35,.65))

png("../Plots/Methods_comparison.png", width=3600, height=2400, res=300)
ggarrange(plotlist=list(p1,p2,p4,
                        p.graph.CLR,p.graph.glasso,
                        p.venn+theme(legend.position="bottom")), 
          labels=c("A","B","C","D","E","F"),
          ncol=3, nrow=2)
dev.off()



# Legend
png(filename="scripts/tmp_files/colorLegend.png", width=1200, height=1200, res=200)
par(mar=c(2,0,2,0))
plot.new()
legend("center",legend=names(colpal),fill=colpal,title="Family",ncol=2)
dev.off()

# Table info
netw.info <- data.frame(
  "Edges"=sapply(list(g.mb, g.gl, g.clr), igraph::ecount),
  
  "Positive"=paste(
    100*sapply(list(g.mb, g.gl, g.clr), 
               function(x) round(sum(igraph::E(x)$weight>0)/igraph::ecount(x),2)),
    "%",sep=""),
  
  "Negative"=paste(
    100*sapply(list(g.mb, g.gl, g.clr), 
               function(x) round(sum(igraph::E(x)$weight<0)/igraph::ecount(x),2)),
    "%",sep="")
)
rownames(netw.info) <- c("MB","GLASSO","CLR")

png(filename="scripts/tmp_files/table.png",width=1200,height=1200,res=300)
ggplot() + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  annotation_custom( gridExtra::tableGrob(netw.info), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
dev.off()

p.mb <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/graph_mb.png") 
p.gl <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/graph_gl.png")
p.clr <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/graph_clr.png") 
p.leg <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/colorLegend.png")
#p.blanck <- ggplot() + theme_void()
p.ven <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/ggvenn.png")
p.inf <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/table.png")

png(filename="outputs/Graph_Comparison.png",width=1600,height=1000,res=300)
gridExtra::grid.arrange(p.mb, p.gl, p.clr, p.ven, p.leg, p.inf, ncol=3)
dev.off()