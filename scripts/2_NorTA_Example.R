library(tidyverse)
library(gridExtra)

# NorTA
rnorm.pos <- mvtnorm::rmvnorm(n=5e2, sigma=matrix(c(1,.8,.8,1),nrow=2))
rnorm.neg <- mvtnorm::rmvnorm(n=5e2, sigma=matrix(c(1,-.8,-.8,1),nrow=2))
rnorm.0 <- mvtnorm::rmvnorm(n=5e2, sigma=matrix(c(1,0,0,1),nrow=2))


unif.pos <- stats::pnorm(rnorm.pos)
unif.neg <- stats::pnorm(rnorm.neg)
unif.0 <- stats::pnorm(rnorm.0)


zinegbin.pos <- apply(unif.pos, 2, function(x)VGAM::qzinegbin(x,munb=30,size=10,pstr0=.25))
htrlnorm.pos <- apply(unif.pos, 2, function(x)ToyModel::qhtrlnorm(x,meanlog=1,b=5,phi=.25))

zinegbin.neg <- apply(unif.neg, 2, function(x)VGAM::qzinegbin(x,munb=30,size=10,pstr0=.25))
htrlnorm.neg <- apply(unif.neg, 2, function(x)ToyModel::qhtrlnorm(x,meanlog=1,b=5,phi=.25))

zinegbin.0 <- apply(unif.0, 2, function(x)VGAM::qzinegbin(x,munb=30,size=10,pstr0=.25))
htrlnorm.0 <- apply(unif.0, 2, function(x)ToyModel::qhtrlnorm(x,meanlog=1,b=5,phi=.25))

p.gen.pos <- ggpubr::ggscatter(as.data.frame(rnorm.pos),
                       x="V1", y="V2",add="reg.line", conf.int=TRUE,
                       add.params = list(color="red",fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), label.x=-3.5, label.y=3, size=6) +
  theme_bw() + ggtitle(label="Normal") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5)) +
  xlim(-3.5,3.5) + ylim(-3.5,3.5)

p.gen.neg <- ggpubr::ggscatter(as.data.frame(rnorm.neg),
                       x="V1", y="V2",add="reg.line", conf.int=TRUE,
                       add.params = list(color="red",fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), label.x=-3.5, label.y=-3, size=6) +
  theme_bw() + ggtitle(label="Normal") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5)) +
  xlim(-3.5,3.5) + ylim(-3.5,3.5)

p.gen.0 <- ggpubr::ggscatter(as.data.frame(rnorm.0),
                       x="V1", y="V2",add="reg.line", conf.int=TRUE,
                       add.params = list(color="red",fill="lightgray")) +
  ggpubr::stat_cor(aes(label=after_stat(r.label)), label.x=-3.5, label.y=-3, size=6) +
  theme_bw() + ggtitle(label="Normal") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5)) +
  xlim(-3.5,3.5) + ylim(-3.5,3.5)

p.unif.pos <- ggplot(as.data.frame(unif.pos), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Normal Quantile") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.unif.neg <- ggplot(as.data.frame(unif.neg), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Normal Quantile") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.unif.0 <- ggplot(as.data.frame(unif.0), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Normal Quantile") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.zinegbin.pos <- ggplot(as.data.frame(zinegbin.pos), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Zero Inflated Negative Binomial") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.zinegbin.neg <- ggplot(as.data.frame(zinegbin.neg), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Zero Inflated Negative Binomial") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.zinegbin.0 <- ggplot(as.data.frame(zinegbin.0), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Zero Inflated Negative Binomial") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.htrlorn.pos <- ggplot(as.data.frame(htrlnorm.pos), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Hurdle Truncated Log-Normal") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.htrlorn.neg <- ggplot(as.data.frame(htrlnorm.neg), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Hurdel Truncated Log-Normal") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

p.htrlorn.0 <- ggplot(as.data.frame(htrlnorm.0), aes(x=V1, y=V2)) +
  geom_point() + theme_bw() + ggtitle("Hurdel Truncated Log-Normal") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust=.5))

png("outputs/NorTA_example.png", width=3600, height=2700, res=300)
grid.arrange(p.gen.pos, p.unif.pos, p.zinegbin.pos, p.htrlorn.pos,
             p.gen.0, p.unif.0, p.zinegbin.0, p.htrlorn.0,
             p.gen.neg, p.unif.neg, p.zinegbin.neg, p.htrlorn.neg,
             nrow=3)
dev.off()