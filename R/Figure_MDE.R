rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))

load(paste0(wd,"Data/Processed/ped_mde.Rdata"))
ncol(all_mde)
mde <- all_mde[3,]

# if a metric includes 0s, add 0.5
includes0<-apply(all_stat,1,function(x)any(x==0))
all_stat2<-all_stat
all_stat2[includes0,]<-all_stat2[includes0,]+0.5

r2s <- apply(all_stat2,1,function(x){
   summary(lm(log(mde)~log(x)))$adj.r.squared
})


setEPS()
pdf(paste0(wd,"Figures/Fig_mde.pdf"), height=6, width=12)
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


  maternities <- as.numeric(all_stat2["maternities",])

  par(mar=c(5,5,3,1))

  plot(mde~maternities,pch=19, ylab="Minimum Detectable Effect",xlab="Number of maternities", xlim=c(0,11000), col= scales::alpha(1,0.5), ylim=c(0,0.2), cex.lab=1.5)

  x<-seq(min(log(maternities)),max(log(maternities)),0.1)
  preds <- as.vector(exp(model.matrix(~ x) %*% coef(lm(log(mde)~log(maternities)))))

  lines(preds~exp(x), col="red")
  mtext("b)", side=3, adj=-0.15,padj=-0.5, cex=2)
  text(8000,0.15, as.expression(bquote(italic(R^2)~"="~.(round(r2s["maternities"],3)))), cex=2)
}
dev.off()


###
# Supplementary Figure with relationships between precision and  all metrics plotted out
###

setEPS()
pdf(paste0(wd,"Figures/FigSM_mde_all.pdf"), height=8, width=14)

{
par(mfrow=c(4,6),mar=c(3.5,3,2,1),mgp=c(2,0.5,0))#mar=c(2,2,3,1))
for(i in 1:23){

  stat <- as.numeric(all_stat2[i,])

 r2 <- summary(lm(log(mde)~log(stat)))$adj.r.squared

 plot(stat,mde,pch=19, col= scales::alpha(1,0.5),xlab=rownames(all_stat)[i], ylab="Minimum Detectable Effect", log="xy", main=paste0("r2=",round(r2,3)))

}}
dev.off()
