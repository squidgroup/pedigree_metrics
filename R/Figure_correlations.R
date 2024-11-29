rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))

all_stat[all_stat==0]<-NA

m <- cor(t(all_stat[1:24,]),use="pairwise.complete.obs")
corrplot(m, order = 'hclust', hclust.method="complete", addrect = 6, diag = FALSE, method = 'square')#, tl.pos = 'd')


# all_stat2<-all_stat

# all_stat2[4:11,] <- all_stat2[4:11,]/((all_stat2[1,]*(all_stat2[1,]-1))/2)
# m2 <- cor(t(all_stat2[1:24,]),use="pairwise.complete.obs")
# corrplot(m2, order = 'hclust', hclust.method="complete", addrect = 6, diag = FALSE, method = 'square', tl.pos = 'd')

## i think the negative relationship between relatedness and pedigree size is maybe expected
# if you add an additional family, then you get some more close relatives, but many more zeros
# also makes sense with more founders 

## relatedness
## you cant get more distant relations without 0.5s
## variance - distribution is 0 bound