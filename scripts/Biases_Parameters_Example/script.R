library(ToyModel)
library(vegan)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

mean_Pielou <- function(x){
  P <- mean( apply(x,1,vegan::diversity) / log(ncol(x)))
  return(round(P,2))
}

cor_D5 <- diag(5); #cor_D10[2,7]<-cor_D10[7,2]<-.8; cor_D10[5,10]<-cor_D10[10,5]<- -.8 
cor_D30 <- diag(30); #cor_D50[21,45]<-cor_D50[45,21]<-.8; cor_D50[6,32]<-cor_D50[32,6]<- -.8
cor_D100 <- diag(100); #cor_D50[21,45]<-cor_D50[45,21]<-.8; cor_D50[6,32]<-cor_D50[32,6]<- -.8

set.seed(10)
toy_D5_P100 <- ToyModel::toy_model(n=10^4, cor=cor_D5, M=1, qdist=qnorm,
                            param=c("mean"=0, "sd"=1), force.positive=T, 
                            method="pearson")

toy_D5_P50 <- ToyModel::toy_model(n=10^4, cor=cor_D5, M=15.5, qdist=qnorm,
                            param=c("mean"=0, "sd"=1), force.positive=T, 
                            method="pearson")


toy_D30_P50 <- ToyModel::toy_model(n=10^4, cor=cor_D30, M=64, qdist=qnorm,
                                     param=c("mean"=0, "sd"=1), force.positive=T, 
                                     method="pearson")

#toy_D100_M100 <- ToyModel::toy_model(n=10^4, cor=cor_D100, M=180, qdist=qnorm,
#                                   param=c("mean"=0, "sd"=.5), force.positive=T, 
#                                   method="pearson")


png(filename="../outputs/biases_example.png", width=8000, height=4400, res=600)

layout(matrix(c(1:12), nrow=4, byrow=T), heights=c(.1,.3,.3,.3))
par(mar=c(0,4,1,4))

plot.new()
#plot.new() rect(xleft=.1, ybottom=.1, xright=.9, ytop=.9, col="lightgray") text(x=.5, y=.5,      paste("Dimension = ",ncol(cor_D10),"\n",            "Magnification = ",1,"\n",            "Pielou Mean = ",mean_Pielou(toy_D10_M1$NorTA)))
text(x=.5, y=.5, "Generated", cex=2, font=2)

plot.new()
#plot.new() rect(xleft=.1, ybottom=.1, xright=.9, ytop=.9, col="lightgray") text(x=.5, y=.5,      paste("Dimension = ",ncol(cor_D10),"\n",            "Magnification = ",1,"\n",            "Pielou Mean = ",mean_Pielou(toy_D10_M1$NorTA)))
text(x=.5, y=.5, "L1", cex=2, font=2)

plot.new()
#plot.new() rect(xleft=.1, ybottom=.1, xright=.9, ytop=.9, col="lightgray") text(x=.5, y=.5,      paste("Dimension = ",ncol(cor_D10),"\n",            "Magnification = ",1,"\n",            "Pielou Mean = ",mean_Pielou(toy_D10_M1$NorTA)))
text(x=.5, y=.5, "CLR", cex=2, font=2)

#plot.new()
#plot.new() rect(xleft=.1, ybottom=.1, xright=.9, ytop=.9, col="lightgray") text(x=.5, y=.5,      paste("Dimension = ",ncol(cor_D10),"\n",            "Magnification = ",1,"\n",            "Pielou Mean = ",mean_Pielou(toy_D10_M1$NorTA)))
# text(x=.5, y=.5,
#      paste("Dimension = ",ncol(cor_D10),"\n",
#            "Magnification = ",1,"\n",
#            "Pielou Mean = ",mean_Pielou(toy_D10_M1$NorTA)))
plot(toy_D5_P100, "Normal", vertex.label=NA)
text(x=-1,y=-1,"A", cex=2)
plot(toy_D5_P100, "L1", vertex.label=NA)
plot(toy_D5_P100, "CLR", vertex.label=NA)

#plot.new()
#plot.new() rect(xleft=.1, ybottom=.1, xright=.9, ytop=.9, col="lightgray") text(x=.5, y=.5,      paste("Dimension = ",ncol(cor_D10),"\n",            "Magnification = ",1,"\n",            "Pielou Mean = ",mean_Pielou(toy_D10_M1$NorTA)))
# text(x=.5, y=.5,
#      paste("Dimension = ",ncol(cor_D10),"\n",
#            "Magnification = ",20,"\n",
#            "Pielou Mean = ",mean_Pielou(toy_D10_M20$NorTA)))
plot(toy_D5_P50, "Normal", vertex.label=NA)
text(x=-1,y=-1,"B", cex=2)
plot(toy_D5_P50, "L1", vertex.label=NA)
plot(toy_D5_P50, "CLR", vertex.label=NA)

# plot.new()
# plot.new() rect(xleft=.1, ybottom=.1, xright=.9, ytop=.9, col="lightgray") text(x=.5, y=.5,      paste("Dimension = ",ncol(cor_D10),"\n",            "Magnification = ",1,"\n",            "Pielou Mean = ",mean_Pielou(toy_D10_M1$NorTA)))
# text(x=.5, y=.5,
#      paste("Dimension = ",ncol(cor_D50),"\n",
#            "Magnification = ",20,"\n",
#            "Pielou Mean = ",mean_Pielou(toy_D50_M20$NorTA)))
plot(toy_D30_P50, "Normal", vertex.label=NA)
text(x=-1,y=-1,"C", cex=2)
plot(toy_D30_P50, "L1", vertex.label=NA)
plot(toy_D30_P50, "CLR", vertex.label=NA)

par(xpd = NA)
abline(h = 1.2, col=rgb(0,0,0,.5))
abline(h = 3.7, col=rgb(0,0,0,.5))
abline(h = 6.3, col=rgb(0,0,0,.5))

dev.off()

mean_Pielou(toy_D5_M1$NorTA)
mean_Pielou(toy_D5_M20$NorTA)
mean_Pielou(toy_D30_M100$NorTA)
