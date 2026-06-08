rm(list=ls())

library(corrplot)

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))

# all_stat[all_stat==0]<-NA

# log all metrics, if a metric includes 0s, add 0.5
includes0<-apply(all_stat,1,function(x)any(x==0))
all_stat2<-all_stat
all_stat2[includes0,]<-log(all_stat2[includes0,]+0.5)


m <- cor(t(all_stat2),use="pairwise.complete.obs")
colnames(m) <- rownames(m) <- gsub("pairwise relatedness","r",rownames(m)) 
colnames(m) <- rownames(m) <- gsub("mean maternal sibship size", "mat. sib. size",rownames(m))
colnames(m) <- rownames(m) <- gsub("mean paternal sibship size", "pat. sib. size",rownames(m))
colnames(m) <- rownames(m) <- gsub("paternal","pat.",rownames(m))
colnames(m) <- rownames(m) <- gsub("maternal","mat.",rownames(m))
colnames(m) <- rownames(m) <- gsub("maximum","max.",rownames(m))

# corrplot(m, order = 'FPC', diag = FALSE, method = 'square')#, tl.pos = 'd')

setEPS()
pdf(paste0(wd,"Figures/Fig_cor.pdf"), height=8, width=8)

# layout(matrix(c(1,2),ncol=2),width=c(1,2))
# plot(NA)

corrplot(m, order = 'hclust', hclust.method="complete", addrect = 8, diag = FALSE, method = 'square',rect.lwd = 4)
	# ,rect.col=palette.colors()[c(2,6,7,3,5,4)],rect.lwd = 4,tl.col=palette.colors()[c(2,6,6,7,7,7,7,7,7,7,7,7,7,7,7,3,3,3,3,3,5,4,4)])#, tl.pos = 'd')
#palette.colors()[c(2,6,7,3,5,4)] c("magenta1")

# lines(c(0.5,0.5),c(0.5,2.5), lwd=4, col=palette.colors()[4])
# lines(c(0.5,0.5),c(2.5,3.5), lwd=4, col=palette.colors()[5])
# lines(c(0.5,0.5),c(3.5,8.5), lwd=4, col=palette.colors()[3])
# lines(c(0.5,0.5),c(8.5,20.5), lwd=4, col=palette.colors()[7])
# lines(c(0.5,0.5),c(20.5,22.5), lwd=4, col=palette.colors()[6])
# lines(c(0.5,0.5),c(22.5,23.5), lwd=4, col=palette.colors()[2])

# lines(c(0.5,1.5),c(23.5,23.5), lwd=4, col=palette.colors()[2])
# lines(c(1.5,3.5),c(23.5,23.5), lwd=4, col=palette.colors()[6])
# lines(c(3.5,15.5),c(23.5,23.5), lwd=4, col=palette.colors()[7])
# lines(c(15.5,20.5),c(23.5,23.5), lwd=4, col=palette.colors()[3])
# lines(c(20.5,21.5),c(23.5,23.5), lwd=4, col=palette.colors()[5])
# lines(c(21.5,23.5),c(23.5,23.5), lwd=4, col=palette.colors()[4])

# lines(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5),c(0.5,2.5,3.5,8.5,20.5,22.5,23.5), lwd=4, col=palette.colors()[c(2,6,7,3,5,4)])


dev.off()

