rm(list=ls())

library(xtable)

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims_low.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims_high.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims.Rdata"))

all_sim_med <- sqrt(all_sim[2,])
all_sim_low <- sqrt(all_sim_low[2,])
all_sim_high <- sqrt(all_sim_high[2,])


includes0<-apply(all_stat,1,function(x)any(x==0))
all_stat2<-all_stat
all_stat2[includes0,]<-all_stat2[includes0,]+0.5

r2s <- t(apply(all_stat2,1,function(x){
  c(
    low=summary(lm(log(all_sim_low)~log(x)))$adj.r.squared,
    med=summary(lm(log(all_sim_med)~log(x)))$adj.r.squared,
    high=summary(lm(log(all_sim_high)~log(x)))$adj.r.squared
  )
}))

xtable(r2s, type = "latex")

setEPS()
pdf(paste0(wd,"Figures/FigSM_h2.pdf"), height=10, width=5)
{
par(mar=c(0,0,0,0), mfrow=c(1,1))
plot(NA, xlim=c(-5,8),ylim=c(0,26))

text(rep(2,23),23:1,rownames(r2s)[1:23], adj=1)
cols <- viridis::viridis(21)[21:1]#hcl.colors(41,palette="TealRose")[41:1]

text(rep(3,23),23:1,round(r2s[,1],3)[1:23], col=cols[as.numeric(cut(r2s[,1],seq(-0.05,1,0.05)))])
text(rep(5,23),23:1,round(r2s[,2],3)[1:23], col=cols[as.numeric(cut(r2s[,2],seq(-0.05,1,0.05)))])
text(rep(7,23),23:1,round(r2s[,3],3)[1:23], col=cols[as.numeric(cut(r2s[,3],seq(-0.05,1,0.05)))])

abline(v=2.25, col="grey")
abline(h=24, col="grey")
text(c(5),c(26),c(expression(italic(h^2))))
text(c(3,5,7),c(25),c("0.1","0.3","0.5"))
text(c(2),c(25),c("Metric"), adj=1)
}

dev.off()