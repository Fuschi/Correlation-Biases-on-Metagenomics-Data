library(tidyverse)
library(SpiecEasi)

# Load custom function
source("R/CLR.R")
source("R/TRIU.R")
source("R/LAYOUT_SIGNED.R")



# READ AND FILTER DATA
#------------------------------------------------------------------------------#

# Read HMP2
otu <- readRDS("data/otu_HMP2.rds")
meta <- readRDS("data/meta_HMP2.rds")
taxa <- readRDS("data/taxonomy.rds")

# Select samples belonging to 69-001 subject in health status
otu.69001.H <- otu[meta$SubjectID=="69-001" & meta$CL4_2=="Healthy" , ]

# Remove rarest OTUs using prevalence and median of non-zero values
otu.filt <- otu.69001.H[, colSums(otu.69001.H>0)/nrow(otu.69001.H) >= .33]
otu.filt <- otu.filt[, apply(otu.filt,2,function(x) median(x[x>0])>=5)]
taxa.filt <- taxa[colnames(otu.filt), ]



# COMPUTES ALL METHODS AND PREPARE RESULTS FOR PLOTS
#------------------------------------------------------------------------------#

# Spiec-Easi MB
res.mb <- SpiecEasi::spiec.easi(data=otu.filt, method='mb',
                                lambda.max=.5, lambda.min.ratio=.5,
                                pulsar.params=list(ncores=6, thresh=0.05))
adj.mb <- as.matrix(SpiecEasi::symBeta(as.matrix(getOptBeta(res.mb)), 
                    mode='maxabs'))
colnames(res.mb$est$data) -> colnames(adj.mb) -> rownames(adj.mb)
adj.mb[adj.mb>0] <- 1
adj.mb[adj.mb<0] <- -1

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
res.mb.clr <- triu(res.clr*abs(adj.mb)); res.mb.clr <- res.mb.clr[res.mb.clr!=0]
res.gl.clr <- triu(res.clr*abs(adj.gl)); res.gl.clr <- res.gl.clr[res.gl.clr!=0]

df.mb <- data.frame("method"=rep("clr",length(triu(res.clr))),"value"=triu(res.clr))
df.mb <- rbind(df.mb, data.frame("method"=rep("mb",length(res.mb.clr)),"value"=res.mb.clr))

df.gl <- data.frame("method"=rep("clr",length(triu(res.clr))),"value"=triu(res.clr))
df.gl <- rbind(df.gl, data.frame("method"=rep("gl",length(res.gl.clr)),"value"=res.gl.clr))



# PLOTS
#------------------------------------------------------------------------------#

# SparCC
p1 <- ggpubr::ggscatter(data.frame("PearsonCLR"=triu(res.clr),
                                   "SparCC"=triu(res.cc)),
                        x="PearsonCLR", y="SparCC",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), 
                   label.x=.45, label.y=-.25, size=6) +
  theme_bw() +
  annotation_custom(grid::textGrob(label="1)", gp=grid::gpar(cex=1.5),
                                   x=unit(0.05,"npc"), y=unit(0.95,"npc"))) +
  ggtitle("SparCC") + xlab("Pearson+CLR") +
  theme(plot.title = element_text(hjust = 0.5))

#Rho
p2 <- ggpubr::ggscatter(data.frame("PearsonCLR"=triu(res.clr),
                                   "Rho"=triu(res.rho)),
                        x="PearsonCLR", y="Rho",
                        add="reg.line", conf.int=TRUE,
                        add.params = list(color="red", fill="lightgray")) +
  ggpubr::stat_cor(aes(label = after_stat(r.label)), label.x=.45, label.y=-.25, size=6) +
  theme_bw() +
  annotation_custom(grid::textGrob(label="2)", gp=grid::gpar(cex=1.5),
                                   x=unit(0.05,"npc"), y=unit(0.95,"npc"))) +
  ggtitle("Rho") + xlab("Pearson+CLR") +
  theme(plot.title = element_text(hjust = 0.5))

#MB
p3.all<-ggplot(df.mb, aes(x=value, fill=method, color=method)) +
  geom_histogram(position="identity", alpha=.5, breaks=seq(-1,1,by=.1)) +
  geom_hline(yintercept=200, linetype="twodash", color="black", linewidth=1) +
  xlim(c(-1,1)) +
  theme_bw() +
  scale_color_manual(values=c("steelblue1", "darkred"))+
  scale_fill_manual(values=c("steelblue1", "darkred"))+
  ggtitle("MB") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")

p3.zoom <- ggplot(df.mb, aes(x=value, fill=method, color=method)) +
  geom_histogram(position="identity", alpha=.5, breaks=seq(-1,1,by=.1)) +
  coord_cartesian(ylim=c(0, 200)) +
  theme_bw() +
  scale_color_manual(values=c("steelblue1", "darkred"))+
  scale_fill_manual(values=c("steelblue1", "darkred")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme_void() +
  theme(legend.position = "none")

p3 <- p3.all +
  annotation_custom(ggplotGrob(p3.zoom), xmin=-1.1, xmax=-.3, 
                    ymin=1400, ymax = 1850) +
  annotation_custom(grid::textGrob(label="3)", gp=grid::gpar(cex=1.5),
                                   x=unit(0.92,"npc"), y=unit(0.95,"npc")))

#GLASSO
p4.all<-ggplot(df.gl, aes(x=value, fill=method, color=method)) +
  geom_histogram(position="identity", alpha=.5, breaks=seq(-1,1,by=.1)) +
  geom_hline(yintercept=200, linetype="twodash", color="black", linewidth=1) +
  xlim(c(-1,1)) +
  theme_bw() +
  scale_color_manual(values=c("steelblue1", "forestgreen"))+
  scale_fill_manual(values=c("steelblue1", "forestgreen")) +
  ggtitle("GLASSO") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")

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
  theme(legend.position="none")

p4 <- p4.all +
  annotation_custom(ggplotGrob(p4.zoom), xmin=-1.1, xmax=-.3, 
                    ymin=1400, ymax = 1850) +
  annotation_custom(grid::textGrob(label="4)", gp=grid::gpar(cex=1.5),
                                   x=unit(0.92,"npc"), y=unit(0.95,"npc")))



# Save output
#------------------------------------------------------------------------------#
png(filename="outputs/Method_Comparison.png",width=2400,height=2400,res=300)
gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()




# GRAPHS COMPARISON
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Networks
g.mb <- igraph::graph_from_adjacency_matrix(adj.mb, mode="undirected",
                                            weighted=TRUE)

g.gl <- igraph::graph_from_adjacency_matrix(adj.gl, mode="undirected",
                                            weighted=TRUE)

g.clr <- igraph::graph_from_adjacency_matrix(adj.clr, mode="undirected",
                                             weighted=TRUE)


# Venn Diagramm GGplot
png("scripts/tmp_files/ggvenn.png",width=1200, height=1200, res=300)
ggVennDiagram::ggVennDiagram(x=list(
  "MB"=paste(igraph::as_edgelist(g.mb)[,1],"-",igraph::as_edgelist(g.mb)[,2],sep=""),
  "GLASSO"=paste(igraph::as_edgelist(g.gl)[,1],"-",igraph::as_edgelist(g.gl)[,2],sep=""),
  "CLR"=paste(igraph::as_edgelist(g.clr)[,1],"-",igraph::as_edgelist(g.clr)[,2],sep="")
), 
label_alpha=0) +
  ggplot2::scale_fill_gradient("Shared \n Links",low="white",high = "red") +
  ggplot2::scale_color_manual(values = c("gray20","gray20","gray20"))
dev.off()

# Plot Networks
#------------------------------------------------------------------------------#

# Generate colormap for taxonomy at family level
set.seed(42)
colpal <- rownames(qualpalr::qualpal(n=length(unique(taxa.filt[,"family"])))$RGB)
names(colpal) <- unique(taxa.filt[,"family"])

# Vertex size proportional to clr abundances
vertex.size <- colMeans(CLR(otu.filt) - min(CLR(otu.filt)))

# MB
png(filename="scripts/tmp_files/graph_mb.png",width=1200, height=1200, res=300)
par(mar=c(0,0,2,0))
set.seed(42)
plot(g.mb, vertex.label=NA, main="MB",
     vertex.color=colpal[taxa.filt[,"family"]],
     vertex.size=vertex.size,
     edge.color=ifelse(igraph::E(g.mb)$weight>0,rgb(0,0,1),rgb(1,0,0)),
     edge.width=.5,
     layout=LAYOUT_SIGNED(g.mb))
dev.off()

# GLASSO
png(filename="scripts/tmp_files/graph_gl.png",width=1200, height=1200, res=300)
par(mar=c(0,0,2,0))
set.seed(42)
plot(g.gl, vertex.label=NA, main="GLASSO",
     vertex.color=colpal[taxa.filt[,"family"]],
     vertex.size=vertex.size,
     edge.color=ifelse(igraph::E(g.gl)$weight>0,rgb(0,0,1),rgb(1,0,0)),
     edge.width=.5,
     layout=LAYOUT_SIGNED(g.gl))
dev.off()

# CLR
png(filename="scripts/tmp_files/graph_clr.png",width=1200, height=1200, res=300)
par(mar=c(0,0,2,0))
set.seed(42)
plot(g.clr, vertex.label=NA, main="Pearson+CLR",
     vertex.color=colpal[taxa.filt[,"family"]],
     vertex.size=vertex.size,
     edge.color=ifelse(igraph::E(g.clr)$weight>0,rgb(0,0,1),rgb(1,0,0)),
     edge.width=.5,
     layout=LAYOUT_SIGNED(g.clr))
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