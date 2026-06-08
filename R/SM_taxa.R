rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))

files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
id <- gsub("\\.|[a-z]","",files)

studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))
dryad_ped_info <- read.csv(paste0(wd,"Data/Raw/dryad_ped_info.csv"))
all_ped_info <- rbind(studies_info[,c("id","species_english","taxonomic.class")], dryad_ped_info[,c("id","species_english","taxonomic.class")])

si <- all_ped_info[match(id,all_ped_info$id),c("id","species_english","taxonomic.class")]
si$n <- as.numeric(all_stat[1,files])

si$taxa <- ifelse(si$taxonomic.class %in% c("bird","mammal"),si$taxonomic.class,"other")

range(si$n)


# n_taxa<- tapply(si$species_english,si$taxa,function(x) length(unique(x)))
# n_taxa2 <- paste0(names(n_taxa)," (", n_taxa,")")
taxa_names <- aggregate(n~species_english+taxa,si,length)
taxa_names2 <- paste0(taxa_names$species_english," (", taxa_names$n,")")


setEPS()
pdf(paste0(wd,"Figures/FigSM_taxa.pdf"), height=12, width=8)

par(mar=c(4,12,1,1))
boxplot(n~paste(taxa,species_english),si, las=1, at=c(42:25, 20:7, 2:1),names=taxa_names2, ylab="", range=0, horizontal = TRUE, xlim=c(1,45), log="x")

boxplot(n~taxa,si,add=TRUE, at= c(44,22,4), col=palette.colors(4)[2:4], xaxt="n", cex.axis=1.5, las=2, range=0, horizontal = TRUE, boxwex=2)
abline(h=c(24,6),lwd=2)

dev.off()
