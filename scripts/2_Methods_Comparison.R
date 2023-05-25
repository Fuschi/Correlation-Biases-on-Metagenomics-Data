library(tidyverse)
library(SpiecEasi)

# Load custom function
source("R/CLR.R")
source("R/TRIU.R")



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
res.rho <- propr(counts=otu.filt, metric="rho")@matrix

# Pearson+CLR
res.clr <- cor(CLR(otu.filt), method="pearson")
diag(res.clr) <- 0
adj.clr.d <- corr.p(r=res.clr, n=nrow(otu.filt),
                    adjust="bonferroni", alpha=.05, ci=FALSE)

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