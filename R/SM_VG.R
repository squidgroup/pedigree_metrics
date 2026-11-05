rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

# save(all_sim_low,file=paste0(wd,"Data/Processed/ped_sims_low.Rdata"))

load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims_low.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims_high.Rdata"))
load(paste0(wd,"Data/Processed/ped_h2_se.Rdata"))
load(paste0(wd,"Data/Processed/ped_h2_se_simple.Rdata"))

all_sim_med <- sqrt(all_sim[2,])
all_sim_low <- sqrt(all_sim_low[2,])
all_sim_high <- sqrt(all_sim_high[2,])

all_se_low<-sqrt(all_sv[,1])
all_se_med<-sqrt(all_sv[,2])
all_se_high<-sqrt(all_sv[,3])

simple_se <- sqrt(all_h2_se_simple)

r2_low <- summary(lm(log(all_sim_low)~log(all_se_low)))$adj.r.squared
r2_med <- summary(lm(log(all_sim_med)~log(all_se_med)))$adj.r.squared
r2_high <- summary(lm(log(all_sim_high)~log(all_se_high)))$adj.r.squared

r2_simple_low <- summary(lm(log(all_sim_low)~log(simple_se)))$adj.r.squared
r2_simple_med <- summary(lm(log(all_sim_med)~log(simple_se)))$adj.r.squared
r2_simple_high <- summary(lm(log(all_sim_high)~log(simple_se)))$adj.r.squared


setEPS()
pdf(paste0(wd,"Figures/FigSM_VG.pdf"), height=12, width=8)
{
par(mfrow=c(3,2), mar=c(5,6,5,1), cex.main=2, cex.lab=1.5)

	plot(simple_se,all_sim_low,pch=19, log="xy", xlab=expression(italic(SE[hat(h^2)] ^{"  "~(1)} )), ylab=expression("Expected SE of"~italic(hat(h^2))),ylim=range(c(sqrt(all_h2_se_simple),sqrt(all_sim))), col= scales::alpha(1,0.5), main=expression(italic(h^2)~"="~0.1))
abline(0,1)
 text(0.07,0.005,as.expression(bquote(italic(R^2)~"="~.(round(r2_simple_low,3)))), cex=1.5)


plot(all_se_low,all_sim_low,log="xy",pch=19, xlab=expression(italic(SE[hat(h^2)] ^{"  "~(2)} )), ylab=expression("Expected SE of"~italic(hat(h^2))), col= scales::alpha(hcl.colors(2,palette="TealRose")[1],0.7), main=expression(italic(h^2)~"="~0.1))
abline(0,1)
 text(0.1,0.01,as.expression(bquote(italic(R^2)~"="~.(round(r2_low,3)))), cex=1.5)



plot(simple_se,all_sim_med,pch=19, log="xy", xlab=expression(italic(SE[hat(h^2)] ^{"  "~(1)} )), ylab=expression("Expected SE of"~italic(hat(h^2))),ylim=range(c(sqrt(all_h2_se_simple),sqrt(all_sim))), col= scales::alpha(1,0.5), main=expression(italic(h^2)~"="~0.3))
abline(0,1)
 text(0.07,0.005,as.expression(bquote(italic(R^2)~"="~.(round(r2_simple_med,3)))), cex=1.5)

plot(all_se_med,all_sim_med,log="xy",pch=19, xlab=expression(italic(SE[hat(h^2)] ^{"  "~(2)} )), ylab=expression("Expected SE of"~italic(hat(h^2))), col= scales::alpha(hcl.colors(2,palette="TealRose")[1],0.7), main=expression(italic(h^2)~"="~0.3))
abline(0,1)
 text(0.1,0.018,as.expression(bquote(italic(R^2)~"="~.(round(r2_med,3)))), cex=1.5)



plot(simple_se,all_sim_high,pch=19, log="xy", xlab=expression(italic(SE[hat(h^2)] ^{"  "~(1)} )), ylab=expression("Expected SE of"~italic(hat(h^2))),ylim=range(c(sqrt(all_h2_se_simple),sqrt(all_sim))), col= scales::alpha(1,0.5), main=expression(italic(h^2)~"="~0.5))
abline(0,1)
 text(0.07,0.005,as.expression(bquote(italic(R^2)~"="~.(round(r2_simple_high,3)))), cex=1.5)

plot(all_se_high,all_sim_high,log="xy",pch=19, xlab=expression(italic(SE[hat(h^2)] ^{"  "~(2)} )), ylab=expression("Expected SE of"~italic(hat(h^2))), col= scales::alpha(hcl.colors(2,palette="TealRose")[1],0.7), main=expression(italic(h^2)~"="~0.5))
abline(0,1)
 text(0.1,0.02,as.expression(bquote(italic(R^2)~"="~.(round(r2_high,3)))), cex=1.5)


}

dev.off()