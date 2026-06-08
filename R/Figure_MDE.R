rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))

load(paste0(wd,"Data/Processed/ped_mde.Rdata"))
all_mde
all_sim_se <- sqrt(all_sim[2,])


# if a metric includes 0s, add 0.5
includes0<-apply(all_stat,1,function(x)any(x==0))
all_stat2<-all_stat
all_stat2[includes0,]<-all_stat2[includes0,]+0.5


r2s <- apply(all_stat2,1,function(x){
   summary(lm(log(all_sim_se)~log(x)))$adj.r.squared
})


setEPS()
pdf(paste0(wd,"Figures/Fig_precision_mat.pdf"), height=5, width=7)
{

maternities <- as.numeric(all_stat2["maternities",])
  par(mfrow=c(1,1), mar=c(5,5,1,1))

plot(maternities,sqrt(all_sim[2,]),pch=19, ylab="Expected SE when h2=0.3",xlab="Number of maternities", xlim=c(0,11000), col= scales::alpha(1,0.5), ylim=c(0,0.2))

x<-seq(min(log(maternities)),max(log(maternities)),0.1)
preds <- as.vector(exp(model.matrix(~ x) %*% coef(lm(log(sqrt(all_sim[2,]))~log(maternities)))))

lines(preds~exp(x), col="red")

}
dev.off()


# plot(all_stat[2,],all_mde[3,],pch=19, ylab="Minimum detectable effect",xlab="Number of maternities", xlim=c(0,11000))