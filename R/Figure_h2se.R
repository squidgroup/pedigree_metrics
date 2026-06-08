rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_h2_se.Rdata"))
load(paste0(wd,"Data/Processed/ped_h2_se_simple.Rdata"))



all_sim_sv <- all_sim[2,]

r2 <- summary(lm(log(sqrt(all_sim_sv))~log(sqrt(all_sv[,2]))))$adj.r.squared
r2_simple <- summary(lm(log(sqrt(all_sim_sv))~log(sqrt(all_h2_se_simple))))$adj.r.squared


setEPS()
pdf(paste0(wd,"Figures/Fig_h2_se.pdf"), height=5, width=10)

{
par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.5)
plot(sqrt(all_h2_se_simple),sqrt(all_sim_sv),pch=19, log="xy", xlab=expression(italic(SE[h^2]^{(1)} )), ylab="Expected SE",ylim=range(c(sqrt(all_h2_se_simple),sqrt(all_sim_sv))), col= scales::alpha(1,0.5))
 text(0.07,0.005,as.expression(bquote(italic(R^2)~"="~.(round(r2_simple,3)))), cex=1.5)
abline(0,1)
  mtext("a)", side=3, adj=-0.2,padj=1, cex=2)

plot(sqrt(all_sv[,2]),sqrt(all_sim_sv),pch=19, log="xy", xlab=expression(italic(SE[h^2]^{(2)} )), ylab="Expected SE", col= scales::alpha(1,0.5))
 text(0.1,0.018,as.expression(bquote(italic(R^2)~"="~.(round(r2,3)))), cex=1.5)
abline(0,1)
  mtext("b)", side=3, adj=-0.2,padj=1, cex=2)

}
dev.off()
