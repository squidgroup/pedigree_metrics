rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

####

load(paste0(wd,"Data/Processed/ped_h2_se.Rdata"))
load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(file=paste0(wd,"Data/Processed/ped_sims.Rdata"))
# load(paste0(wd,"/Data/Processed/GRM_lambdas.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_dat.Rdata"))

load(file=paste0(wd,"/Data/Processed/grm_sims.Rdata"))

grm_dat2 <-cbind(as.data.frame(grm_dat),
	col=c("red","blue","green","orange","purple","pink"),
	pch=c(21:25,21)
	)

grm_dat3<- cbind(grm_dat2[rownames(grm_sim),],grm_sim)


all_se <- sqrt(all_sv[,2])
n<-as.numeric(all_stat[1,])
all_se_sim <- sqrt(all_sim[2,])

setEPS()
pdf(paste0(wd,"Figures/Fig_grm.pdf"), height=5, width=10)
{
	par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.5)

	plot(all_se_sim~n, log="xy", ylim=range(c(all_se,grm_dat2[,"se"])), col=scales::alpha(1,0.5), pch=19, ylab="Expected Standard Error", xlab="N")
	points(sqrt(var)~n,grm_dat3, pch=grm_dat3[,"pch"], bg=grm_dat3[,"col"], col=grm_dat3[,"col"], cex=1.5)

	legend(
		"topright",
		c("pedigrees",rownames(grm_dat)),
		pch=c(19,grm_dat2[,"pch"]),
		pt.bg=c(scales::alpha(1,0.5),grm_dat2[,"col"]),
		col=c(scales::alpha(1,0.5),grm_dat2[,"col"])
	)

  mtext("a)", side=3, adj=-0.2,padj=1, cex=2)




	plot(all_se_sim~all_se, log="xy", ylim=range(c(all_se,grm_dat[,"se"])),, xlim=range(c(all_se,grm_dat[,"se"])), col=scales::alpha(1,0.5), pch=19, ylab="Expected Standard Error", xlab=expression(italic(SE[h^2]^{(2)} )))
	
	points(sqrt(var)~se,grm_dat3, pch=grm_dat3[,"pch"], bg=grm_dat3[,"col"], col=grm_dat3[,"col"], cex=1.5)
	abline(0,1)

  mtext("b)", side=3, adj=-0.2,padj=1, cex=2)


}
dev.off()

# grm_vars<-c(var(dragon_grm[lower.tri(dragon_grm)]),
# var(koala_grm[lower.tri(koala_grm)]),
# var(taz_devil[lower.tri(taz_devil)]),
# var(s_grm[lower.tri(s_grm)])
# )
# hist(c(all_stat2[24,]))
# abline(v=grm_vars)

