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

si$taxa <- ifelse(si$taxonomic.class %in% c("bird","mammal"),si$taxonomic.class,"other")

range(si$n)

par(mar=c(4,4,1,1))
boxplot(n~taxonomic.class,si)

n_taxa<- tapply(si$species_english,si$taxa,function(x) length(unique(x)))
taxa_names <- aggregate(n~species_english+taxa,si,length)$species_english


par(mar=c(12,4,1,1))
boxplot(n~paste(taxa,species_english),si, las=2, at=c(2: 19, 22:35, 38:39),names=taxa_names)

boxplot(n~taxa,si,add=TRUE, at= c(1,21,37), col=2:5, yaxt="n")
abline(v=c(20,37))



table(paste(si$taxonomic.class,si$species_english))


par(mar=c(12,4,1,1))
plot(table(si$species), las=2, ylab="Number of pedigrees")


