# Load packages #
library(dplyr)
library(plyr)
require(data.table)

# Read Files #
mapped_read <- read.table("All.sample.mapped_read", row.names=1)
chm13 <- sum(read.table("/gmi-l1/_90.User_Data/juhyunk/Files/STAR/CHM13_GRCh38/CHM13/chrLength.txt")[,1])
genome_len_db <- read.csv("/gmi-l1/_90.User_Data/juhyunk/Files/BlastDB/virus_only_ref/virus.virome.DB/virus_genome_for_first_blast.db.gene.length", header=F, sep='~')
virusCluster <- read.csv("/gmi-l1/_90.User_Data/juhyunk/codes/virome_RNA-seq/virus.Cluster.txt", header=T, sep="|")

# Modify Inputs #
colnames(genome_len_db) <- c("taxonomy_id","len")
genome_len_db$taxonomy_id <- sub(".*\\|","",sub("_","|",sub("^AC_","AC",sub("^NC_","NC",genome_len_db$taxonomy_id))))
genome_len_db<-arrange(genome_len_db, taxonomy_id)

virusCluster$virusName <- gsub(" ","_",sub(".*\\|","",sub(" ","|",virusCluster$virusName)))

taxon_name <- system("cat */BLAST/*.second_blast.firstHit.virusOnly.pairedCountBinningFiltered.txt| cut -f 12 | sort | uniq | sed -e 's/ /_/g'", intern=T)
taxon_name <- c(taxon_name, virusCluster[virusCluster$clusterName %in% taxon_name[startsWith(taxon_name,"cluster")],]$virusName)
genome_len_db <- genome_len_db[genome_len_db$taxonomy_id %in% taxon_name,]



for (virusCluster_num in 1:dim(virusCluster)[1]) {
		genome_len_db$taxonomy_id <- sub(virusCluster[virusCluster_num,1],virusCluster[virusCluster_num,2],genome_len_db$taxonomy_id)
}
genome_len_db <- as.data.table(genome_len_db)
genome_len_db <- as.data.frame(genome_len_db[genome_len_db[, .I[len == max(len)], by=taxonomy_id]$V1])

genome_len <- as.matrix(1/genome_len_db$len)

for (virusCluster_num in 1:dim(virusCluster)[1]) {
	taxon_name <- unique(sub(virusCluster[virusCluster_num,1],virusCluster[virusCluster_num,2],taxon_name))
}

taxon_name <- as.data.frame(taxon_name)
colnames(taxon_name) <- "taxonomy_id"
taxon_name <- arrange(taxon_name, taxonomy_id)

norm_table <- data.frame(matrix(ncol=dim(mapped_read)[1], nrow=dim(taxon_name)[1]))
colnames(norm_table) <- rownames(mapped_read)
rownames(norm_table) <- taxon_name$taxonomy_id



# raw count data 만들기 + human mapped + chm13 genome size normalization
for ( sample_num in 1:dim(mapped_read)[1] ) {
	sample <- rownames(mapped_read)[sample_num]
	#print(sample)
	blast <- as.numeric(system(paste0("ls -l ",sample,"/BLAST/",sample,".second_blast.firstHit.virusOnly.pairedCountBinningFiltered.txt","| cut -d' ' -f 5"),intern=T))
	if (blast > 0 ) {
		blast <- read.table(paste0(sample,"/BLAST/",sample,".second_blast.firstHit.virusOnly.pairedCountBinningFiltered.txt"), header=T, sep='\t')[,12]
	}
	

	if( is.numeric(blast)	) {
		print(sample)
		next
	}

	blast <- gsub(" ", "_",blast)

	blast<-as.data.frame(blast)
	colnames(blast) <- "taxonomy_id"
	blast <- arrange(blast, taxonomy_id)
	blast <- count(blast)

	blast <- t(as.matrix(merge(taxon_name, blast, by="taxonomy_id", all.x=TRUE)[,2]*chm13/as.numeric(mapped_read[sample_num,1])*2))
	blast[is.na(blast)] <- 0
	for ( virus_num in 1:length(genome_len) ) {
	blast[virus_num] <- blast[virus_num]*genome_len[virus_num]
	}
	blast[is.na(blast)] <- 0
	norm_table[,sample_num] <- blast[1,]
}

# 모두 NA인 행 제거
norm_table <- norm_table[!apply(norm_table, 1, function(x) all(x == 0)), ]
norm_table[,"taxonomy_id"] <- rownames(norm_table)

rownames(norm_table) <- NULL
write.table(norm_table,"BLAST.results.YM.txt",sep='\t', quote=FALSE, row.names = F)
