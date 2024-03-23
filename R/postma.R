dd<-read.csv("/Users/joelpick/github/maternal_effects/Data/Raw/postma2014_SM.csv")

names(dd)

table(dd$method.2,dd$publication.year)