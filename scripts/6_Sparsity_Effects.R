library(tidyverse)
library(gridExtra)
library(reshape2)
library(RColorBrewer)


# Set Parameters, CORv=Correlations and PHIv=zero percentage
PHIv <- seq(0,.9,by=.05)
CORv <- seq(-.9,.9,by=.1)


# Logical variable permits to choose if redo all the counts or to select files from cache
# folder
CACHE <- TRUE
stopifnot(is.logical(CACHE))

if(!CACHE){
  
  diffL1 <- matrix(NA,nrow=length(CORv),ncol=length(PHIv),
                 dimnames=list(paste("c",CORv,sep=""),paste("p",PHIv,sep="")))
  diffCLR <- matrix(NA,nrow=length(CORv),ncol=length(PHIv),
                   dimnames=list(paste("c",CORv,sep=""),paste("p",PHIv,sep="")))


  for(PHI in 1:length(PHIv)){
    print(PHI)
    for(COR in 1:length(CORv)){

      corM <- diag(100)
      corM[4,5] <- corM[5,4] <- CORv[COR]
      param <- stats::setNames(c(1,PHIv[PHI]), c("meanlog","phi"))
      toy <- ToyModel::toy_model(n=1e4, cor=corM, M=1, 
                                 qdist=ToyModel::qhtrlnorm,
                                 param=param)

      cor0 <- abs(toy$cor_normal[4,5])
      corL1 <- abs(toy$cor_L1[4,5])
      corCLR <- abs(toy$cor_CLR[4,5])

      diffL1[COR,PHI] <- abs(cor0 - corL1)
      diffCLR[COR,PHI] <- abs(cor0 - corCLR)

    }
  }
} else {
  diffL1  <- readRDS("cache/6_Sparsity_Effects/diffL1.rds")
  diffCLR <- readRDS("cache/6_Sparsity_Effects/diffCLR.rds")
}

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

df.L1 <- melt(log10(diffL1)); df.L1$value[df.L1$value< -2] <- -2
p1 <- ggplot(df.L1, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradientn("Absolute\nError",limits=c(-2,0), colours=myPalette(11),
                       breaks=c(-0.1,-1,-1.9), labels=c("1","0.1","< 0.01")) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_x_discrete(breaks=c("c-0.8","c-0.6","c-0.4","c-0.2","c0","c0.2","c0.4","c0.6","c0.8"),
                   labels=c("-0.8","-0.6","-0.4","-0.2","0","0.2","0.4","0.6","0.8")) +
  scale_y_discrete(breaks=c("p0","p0.1","p0.2","p0.3","p0.4","p0.5","p0.6","p0.7","p0.8","p0.9"),
                   labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Zero %") + xlab("Correlation") + 
  ggtitle("L1 Normalization") + theme(plot.title = element_text(hjust = 0.5)) 

df.CLR <- melt(log10(diffCLR)); df.CLR$value[df.CLR$value< -2] <- -2
#png(filename="L1_effects.png",width=1200,height=1200,res=300)
p2 <- ggplot(df.CLR, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradientn("Absolute\nError",limits=c(-2,0), colours=myPalette(11),
                       breaks=c(-0.1,-1,-1.9), labels=c("1","0.1","< 0.01")) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_x_discrete(breaks=c("c-0.8","c-0.6","c-0.4","c-0.2","c0","c0.2","c0.4","c0.6","c0.8"),
                   labels=c("-0.8","-0.6","-0.4","-0.2","0","0.2","0.4","0.6","0.8")) +
  scale_y_discrete(breaks=c("p0","p0.1","p0.2","p0.3","p0.4","p0.5","p0.6","p0.7","p0.8","p0.9"),
                   labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Zero %") + xlab("Correlation") + 
  ggtitle("CLR Normalization") + theme(plot.title = element_text(hjust = 0.5)) 


png(filename="outputs/Sparsity_Effects.png",width=2400, height=1200, res=300)
p <- list(p1,p2)
grid.arrange(  grobs = p,
               layout_matrix = matrix(c(1, 2),ncol=2))
dev.off()