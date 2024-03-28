
wd <- "/Users/joelpick/github/pedigree_metrics/"
##subset to get just journal articles

kruuk <- subset(read.csv(paste0(wd,"Data/Raw/cited_kruuk2004.csv")),Publication.Type=="J")

wilson <- subset(read.csv(paste0(wd,"Data/Raw/cited_wilson2010.csv")),Publication.Type=="J")

yp <- read.csv(paste0(wd,"Data/Raw/young_postma_2023_data.csv"))
am_papers<- aggregate(estimate.ID~study.ID+	year+	reference+	species.english+	species.latin+	taxonomic.class+	population+	country,subset(yp,yp$method.2 %in% c("animal model (MCMC)","animal model (REML)")),length)

names(am_papers)

sum(kruuk$Article.Title %in% wilson$Article.Title)

both<-rbind(kruuk, wilson[!wilson$Article.Title %in%kruuk$Article.Title,])
nrow(both)

nrow(am_papers)

dup<-do.call(rbind,lapply(1:nrow(both), function(i){
	dd <-grep(tolower(both$Article.Title)[i],tolower(am_papers$ref))
	if(length(dd)>0) data.frame(new_list=i,yp=dd)
}))

## these all match
# rbind(both[dup$new_list,c("Article.Title")],am_papers$ref[dup$yp])
# both[dup,c("Authors","Article.Title","Publication.Year")]

to_screen <- both[-(unique(dup$new_list)),c("Authors","Publication.Year","Article.Title","Source.Title","Abstract")]
nrow(to_screen)

names(to_screen) <- c("Author","Publication.Year","Title","Publication.Title","Abstract")

write.csv(to_screen, file=paste0(wd,"Data/Processed/new_to_screen.csv"),row.names=FALSE)

