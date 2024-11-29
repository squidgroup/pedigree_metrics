rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))

files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
id <- gsub("\\.|[a-z]","",files)

studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))
dryad_ped_info <- read.csv(paste0(wd,"Data/Raw/dryad_ped_info.csv"))
all_ped_info <- rbind(studies_info[,c("id","species_english","taxonomic.class")], dryad_ped_info[,c("id","species_english","taxonomic.class")])

si <- all_ped_info[match(id,all_ped_info$id),c("id","species_english","taxonomic.class")]
si$n <- all_stat[1,files]

par(mar=c(4,4,1,1))
boxplot(n~taxonomic.class,si)

par(mar=c(12,4,1,1))
boxplot(n~species_english,si, las=2)

par(mar=c(12,4,1,1))
plot(table(si$species), las=2, ylab="Number of pedigrees")


