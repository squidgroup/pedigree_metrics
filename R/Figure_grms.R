rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

####

load(paste0(wd,"Data/Processed/ped_h2_se.Rdata"))
load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(file=paste0(wd,"Data/Processed/ped_sims.Rdata"))
# load(paste0(wd,"/Data/Processed/GRM_lambdas.Rdata"))

load(file=paste0(wd,"/Data/Processed/grm_roe_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_soay_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_koala_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_dragon_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_td_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_stickleback_dat.Rdata"))



grm_dat <- rbind(
	roe= roe_dat,
	soays = soays_dat,
	koala = koala_dat,
	dragon = dragon_dat,
	taz_devil = taz_devil_dat,
	stickleback = stickleback_dat

	)


####


grm_dat2 <-cbind(as.data.frame(grm_dat),
	col=c("red","blue","green","orange","purple","pink"),
	pch=c(21:25,21)
	)

all_se <- sqrt(all_sv[,2])
n<-as.numeric(all_stat[1,])
all_se_sim <- sqrt(all_sim[2,])

cbind(all_se[n<600&n>400],n[n<600&n>400])


setEPS()
pdf(paste0(wd,"Figures/Fig_grm.pdf"), height=5, width=10)
{
	par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.5)

	plot(all_se_sim~n, log="xy", ylim=range(c(all_se,grm_dat2[,"se"])), col=scales::alpha(1,0.5), pch=19, ylab="Expected Standard Error", xlab="N")
	points(sqrt(var)~n,grm_dat2, pch=grm_dat2[,"pch"], bg=grm_dat2[,"col"], col=grm_dat2[,"col"], cex=1.5)

	legend(
		"topright",
		c("pedigrees",rownames(grm_dat)),
		pch=c(19,grm_dat2[,"pch"]),
		pt.bg=c(scales::alpha(1,0.5),grm_dat2[,"col"]),
		col=c(scales::alpha(1,0.5),grm_dat2[,"col"])
	)

  mtext("a)", side=3, adj=-0.25,padj=1, cex=2)




	plot(all_se_sim~all_se, log="xy", ylim=range(c(all_se,grm_dat[,"se"])),, xlim=range(c(all_se,grm_dat[,"se"])), col=scales::alpha(1,0.5), pch=19, ylab="Expected Standard Error", xlab=expression(italic(SE[h^2]^{(2)} )))
	
	points(sqrt(var)~se,grm_dat2, pch=grm_dat2[,"pch"], bg=grm_dat2[,"col"], col=grm_dat2[,"col"], cex=1.5)
	abline(0,1)

  mtext("b)", side=3, adj=-0.25,padj=1, cex=2)


}
dev.off()

