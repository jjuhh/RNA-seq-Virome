library(ggplot2)
library(data.table)
library(dplyr)

args <- commandArgs()
similarity <- args[1]
tab <- as.data.frame(fread(similarity, sep='|'))

tab <- tab[tab$V1!=tab$V2,] %>% arrange(-V3)
tab[,"bin"]<- log(c(1:dim(tab)[1]))

png("similairity.png")
ggplot(data=tab, aes(x=bin, y=V3)) + geom_point()
dev.off()

