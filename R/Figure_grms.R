rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))


####

load(paste0(wd,"Data/Processed/ped_h2_se.Rdata"))
load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(paste0(wd,"/Data/Processed/GRM_lambdas.Rdata"))

all_se <- sqrt(all_sv[,2])
n<-as.numeric(all_stat[1,])

setEPS()
pdf(paste0(wd,"Figures/Fig_grm.pdf"), height=5, width=6)
{
par(mar=c(4,4,1,1))
plot(all_se~n, log="xy", ylim=range(c(all_se,s_se)), col=scales::alpha(1,0.5), pch=19, ylab="Approximate Standard Error", xlab="N")
points(se~n,grm_dat,pch=c(21:25,21),
	bg=c("red","blue","green","orange","purple","pink"),
	col=c("red","blue","green","orange","purple","pink"),cex=1.5)

legend(
	"topright",
	c("pedigrees",rownames(grm_dat)),#,"Taz. Devil","Stickleback","Dragon","Koala","Roe Deer","Soay Sheep"
	pch=c(19,21:25,21),
	pt.bg=c(scales::alpha(1,0.5),"red","blue","green","orange","purple","pink"),
	col=c(scales::alpha(1,0.5),"red","blue","green","orange","purple","pink")

	)

}
dev.off()

# grm_vars<-c(var(dragon_grm[lower.tri(dragon_grm)]),
# var(koala_grm[lower.tri(koala_grm)]),
# var(taz_devil[lower.tri(taz_devil)]),
# var(s_grm[lower.tri(s_grm)])
# )
# hist(c(all_stat2[24,]))
# abline(v=grm_vars)

