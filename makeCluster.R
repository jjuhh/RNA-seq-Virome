# similarity가 90이상인 것들 뽑아내기
# makeCluster.R

library(ggplot2)
library(dplyr)
library(data.table)

args <- commandArgs()
stretcherResults <- args[1]

tab <- as.data.frame(fread(stretcherResults, sep='|'))

globalSimilarity <- 90

sub <- tab[tab$V3>globalSimilarity,]
sub<- sub[sub$V1!=sub$V2,]

clusterDb <- data.frame()
clusterCount <- 0

for ( taxon1_num in 1:length(unique(sub$V1)) ){ 
	print(unique(sub$V1)[taxon1_num])
	cluster <- c()
	cluster <- c(cluster ,unique(sub$V1)[taxon1_num] )
	cluster <- c(cluster, sub[sub$V1==unique(sub$V1)[taxon1_num] & sub$V3 > globalSimilarity,]$V2)
	if (length(cluster) > 1 ) { 
			clusterCount <- clusterCount+1	
			print(clusterCount)
		for ( taxon2_num in 1:length(cluster) ) {
				print(taxon2_num)
				cluster <- c(cluster, sub[sub$V1==cluster[taxon2_num] & sub$V3 > globalSimilarity,]$V2)
				cluster <- unique(cluster)
		}
	clusterDb <- rbind(clusterDb, data.frame(virusName=cluster, clusterName=rep(paste0("clusterNum_",clusterCount), length(cluster))))
	}
}

duplicate <- rownames(table(clusterDb[,"virusName"])>1)[table(clusterDb[,"virusName"])>1]

cluster<-c()
for ( taxon2_num in 1:length(duplicate) ) {
	print(taxon2_num)
	cluster <- c(cluster, clusterDb[clusterDb$virusName == duplicate[taxon2_num],]$clusterName)
	cluster <- unique(cluster)
}

clusterDb[,"mergedClusterName"] <- clusterDb$clusterName

clusterList <- list()
for (i in 1:length(duplicate) ){
print(duplicate[i])
clusterList[i] <- list(clusterDb[clusterDb==duplicate[i],"clusterName"])
}
clusterList <- unique(clusterList)

clusterList2 <- list()
for (i in 1:length(table(unlist(clusterList))) ) {
print(i)
originalLen <- length(clusterList)
clusterName <- names(sort(table(unlist(clusterList)),decreasing=T)[i])
clusterNameList <- unique(unlist(clusterList[grep(clusterName, clusterList)]))
clusterList2[i] <- list(unique(unlist(clusterList[grep(clusterName, clusterList)])))
for (j in clusterNameList) {
clusterList2[[i]]<- c(clusterList2[[i]],unique(unlist(clusterList[grep(j, clusterList)])))
}
clusterList2[[i]] <- unique(clusterList2[[i]])
}
clusterList2 <- unique(unique(clusterList2))

#### ---------------------------------------------------------------------

clusterList2 <- unique(clusterList2)
originalLen <- length(unique(clusterList))

while (originalLen > length(unique(clusterList2))) {
clusterList <- unique(clusterList2)
originalLen <- length(unique(clusterList))

for (i in 1:length(unique(unlist(clusterList2))) ) {
print(i)
clusterName <- names(sort(table(unlist(clusterList)),decreasing=T)[i])
clusterNameList <- unique(unlist(clusterList[grep(clusterName, clusterList)]))
clusterList2[i] <- list(unique(unlist(clusterList[grep(clusterName, clusterList)])))
for (j in clusterNameList) {
clusterList2[[i]]<- c(clusterList2[[i]],unique(unlist(clusterList[grep(j, clusterList)])))
}
clusterList2[[i]] <- unique(clusterList2[[i]])
}

clusterList2 <- unique(clusterList2)
}

### ----------------------------------------------------------------------
for (i in 1:length(clusterList2)) {
print(i)
clusterDb[clusterDb$clusterName %in% unlist(clusterList2[i]),"mergedClusterName"] <- paste0("mergedCluster_",i)
}

### ---------------------------------------------------------------------


write.table(clusterDb, "virus.Cluster.txt",quote=F, col.names=T, row.names=F, sep="|")
