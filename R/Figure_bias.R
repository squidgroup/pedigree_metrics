rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(paste0(wd,"Data/Processed/ped_mbias.Rdata"))

all_bias <- all_mbias_sims[1,]

includes0<-apply(all_stat,1,function(x)any(x==0))
all_stat2<-all_stat
all_stat2[includes0,]<-all_stat2[includes0,]+0.5


r2s <- apply(all_stat2,1,function(x){
   summary(lm(log(all_bias)~log(x)))$adj.r.squared
})


###############
## bias Figure
###############

setEPS()
pdf(paste0(wd,"Figures/Fig_bias.pdf"), height=6, width=12)
{

  layout(matrix(c(1,2),ncol=2),width=c(1,2))
  par(mar=c(1,1,1,1))
  plot(NA, xlim=c(-3,4),ylim=c(1,25),bty="n", xaxt="n", yaxt="n")

  text(rep(2,23),23:1,names(r2s)[1:23], adj=1)
  text(rep(3,23),23:1,round(r2s,3)[1:23], col=viridis::viridis(21)[21:1][as.numeric(cut(r2s,seq(-0.05,1,0.05)))])
  abline(v=2.25, col="grey")
  abline(h=24, col="grey")
  text(c(3),c(25),c(expression(italic(R^2))), cex=1.2)
  text(c(2),c(25),c("Metric"), adj=1, cex=1.2)
  mtext("a)", side=3, adj=0.,padj=1, cex=2)


  mat_sib <- as.numeric(all_stat2["mean maternal sibship size",])
  par(mar=c(5,5,3,1))

  plot(mat_sib,all_bias,pch=19, ylab=expression(Bias~"in"~italic(h^2)),xlab="Mean maternal sibship size", col= scales::alpha(1,0.5), cex.lab=1.5)
  x<-seq(min(log(mat_sib)),max(log(mat_sib)),0.1)
  preds <- as.vector(exp(model.matrix(~ x) %*% coef(lm(log(all_bias)~log(mat_sib)))))

  lines(preds~exp(x), col="red")
  mtext("b)", side=3, adj=-0.15,padj=-0.5, cex=2)
  text(12,0.15, as.expression(bquote(italic(R^2)~"="~.(round(r2s["mean maternal sibship size"],3)))), cex=2)
}
dev.off()



setEPS()
pdf(paste0(wd,"Figures/FigSM_bias_all.pdf"), height=8, width=14)
{

par(mfrow=c(4,6),mar=c(3.5,3,2,1),mgp=c(2,0.5,0))#mar=c(2,2,3,1))
for(i in 1:23){

  stat <- as.numeric(all_stat2[i,])

  r2 <- summary(lm(log(all_bias)~log(stat)))$adj.r.squared

  plot(stat,all_bias,pch=19, col= scales::alpha(1,0.5),xlab=rownames(all_stat)[i], ylab=expression("Bias in"~italic(h^2)), log="xy", main=as.expression(bquote(italic(R^2)~"="~.(round(r2,3)))))
  }


# r2 <- summary(lm(log(all_bias)~log(sqrt(all_sv))))$adj.r.squared

# # plot(all_sim[2,],all_m,pch=19, log="xy", main="Visscher & Goddard 2015")
# plot(sqrt(all_sv),all_bias,pch=19, log="xy", xlab="V-G SE", ylab="Bias in Va", col= scales::alpha(1,0.5), main=paste0("r2=",round(r2,3)))

}
dev.off()
