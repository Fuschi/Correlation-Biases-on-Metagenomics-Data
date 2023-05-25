library(tidyverse)

# Set Parameters, D=Dimensionality, M=Magnification Factor
D <- seq(5,200,by=5)
M <- c(1:10,seq(12,98,by=2),seq(100,1000,by=50),seq(1500,10000,by=500))

# load TOY_MODEL function
source("R/TOY_MODEL.R")

# Logical variable permits to choose if redo all the counts or to select files from cache
# folder
CACHE <- TRUE
stopifnot(is.logical(CACHE))

# Elaborates Transformations Effects
#------------------------------------------------------------------------------#

if(!CACHE){
  
  effL1.mean  <- matrix(NA_real_,nrow=length(D),ncol=length(M),dimnames=list(paste("d",D,sep=""),paste("m",M,sep="")))
  effCLR.mean <- matrix(NA_real_,nrow=length(D),ncol=length(M),dimnames=list(paste("d",D,sep=""),paste("m",M,sep="")))
  shannon <- matrix(NA_real_,nrow=length(D),ncol=length(M),dimnames=list(paste("d",D,sep=""),paste("m",M,sep="")))

  for(d in 1:length(D)){
    
    print(d)
    for(m in 1:length(M)){
      
      toy <- TOY_MODEL(n=10^4,cor=diag(D[d]),M=M[m],dist="norm",
                       method="pearson",force.positive=TRUE)
      
      tmpL1 <- abs(toy$cor_NorTA - toy$cor_L1)
      tmpCLR <- abs(toy$cor_NorTA - toy$cor_CLR)
      
      effL1.mean[d,m]   <- mean(tmpL1)
      effCLR.mean[d,m]  <- mean(tmpCLR)
      
      shannon[d,m] <- mean(apply(toy$NorTA,1,vegan::diversity))
      
    }
  }
  shannon.norm <- apply(shannon,2,function(x)x/log(D))
  
} else {
  
  effL1.mean <- readRDS("cache/effL1_mean_0.rds")
  effCLR.mean <- readRDS("cache/effCLR_mean_0.rds")
  shannon.norm <- readRDS("cache/shannon_norm.rds")
  
}

# Set Entropy Values
#------------------------------------------------------------------------------#
E <- rev(seq(.025,1,by=.025))

# Order errors respect the mean normailzed entropy
#------------------------------------------------------------------------------#
effL1.final.mean  <- matrix(NA_real_,nrow=length(E),ncol=length(D),
                            dimnames=list(paste("d",D,sep=""),paste("e",E,sep="")))
effCLR.final.mean <- matrix(NA_real_,nrow=length(E),ncol=length(D),
                            dimnames=list(paste("d",D,sep=""),paste("e",E,sep="")))

for(d in 1:length(D)){
  for(e in 1:length(E)){
    
    ei <- shannon.norm[d,]
    idxe <- which(abs(E[e]-ei)==min(abs(E[e]-ei)),arr.ind=T)
    
    effL1.final.mean[d,e] <- effL1.mean[d,idxe]
    effCLR.final.mean[d,e] <- effCLR.mean[d,idxe]
  }
}


# PLOT
#------------------------------------------------------------------------------#
myPalette <- RColorBrewer::brewer.pal(11, "Spectral") %>%
  rev() %>% grDevices::colorRampPalette()
dfL1.mean  <- reshape2::melt(log10(effL1.final.mean))
dfCLR.mean <- reshape2::melt(log10(effCLR.final.mean))

# Threshold at 0.01
dfL1.mean$value[dfL1.mean$value < -2] <- -2
dfCLR.mean$value[dfCLR.mean$value < -2] <- -2

# L1 Effects
#png(filename="outputs/L1_effects.png",width=1200,height=1200,res=300)
pL1 <- ggplot(dfL1.mean, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradientn("Mean Error",limits=c(-2,0), colours=myPalette(11),
                       breaks=c(-0.1,-1,-1.9), labels=c("1","0.1","< 0.01")) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_discrete(breaks=c("e0.025","e0.1","e0.2","e0.3","e0.4","e0.5","e0.6","e0.7","e0.8","e0.9","e1"),
                   labels=c("0.025","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")) +
  scale_x_discrete(breaks=c("d10","d20","d30","d40","d50","d60","d70","d80","d90",
                            "d100","d110","d120","d130","d140","d150","d160","d170","d180","d190","d200"),
                   labels=c("10","20","30","40","50","60","70","80","90",
                            "100","110","120","130","140","150","160","170","180","190","200")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(expression(paste(alpha,"-Heterogeneity"))) + xlab("Dimensions") + 
  ggtitle("L1 Normalization") + theme(plot.title = element_text(hjust = 0.5)) 
#dev.off()

# CLR Effects
#png(filename="CLR_effects.png",width=1200,height=1200,res=300)
pCLR <- ggplot(dfCLR.mean, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradientn("Mean Error",limits=c(-2,0), colours=myPalette(11),
                       breaks=c(-0.1,-1,-1.9), labels=c("1","0.1","< 0.01")) +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size=10)) +
  scale_y_discrete(breaks=c("e0.025","e0.1","e0.2","e0.3","e0.4","e0.5","e0.6","e0.7","e0.8","e0.9","e1"),
                   labels=c("0.025","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")) +
  scale_x_discrete(breaks=c("d10","d20","d30","d40","d50","d60","d70","d80","d90",
                            "d100","d110","d120","d130","d140","d150","d160","d170","d180","d190","d200"),
                   labels=c("10","20","30","40","50","60","70","80","90",
                            "100","110","120","130","140","150","160","170","180","190","200")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(expression(paste(alpha,"-Heterogeneity"))) + xlab("Dimensions") + 
  ggtitle("CLR Normalization") + 
  theme(plot.title = element_text(hjust = 0.5))  + theme(plot.title = element_text(hjust = 0.5))
#dev.off


png(filename="outputs/Normalization_Bias.png",width=2400, height=1200, res=300)
p <- list(pL1,pCLR)
#X11(width=12,height=6)
gridExtra::grid.arrange(  grobs = p,
               layout_matrix = matrix(c(1,2),ncol=2))
dev.off()
