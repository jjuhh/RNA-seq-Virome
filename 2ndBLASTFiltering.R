## second blast filtering ##
# Rscript: 

rm(list=ls())

# Parameters #
args <- commandArgs()
binningCutoff <- args[1]
countCutoff <- args[2]

# Load packages #
library(dplyr)
library(plyr)
require(data.table)

# Read Files #
mapped_read <- read.table("All.sample.mapped_read", row.names=1)
genome_len_db <- read.csv("/gmi-l1/_90.User_Data/juhyunk/Files/BlastDB/virus_only_ref/virus.virome.DB/virus_genome_for_first_blast.db.gene.length", header=F, sep='~', col.names=c("taxonomy_id","len"))
virusCluster <- read.csv("/gmi-l1/_90.User_Data/juhyunk/codes/virome_RNA-seq/virus.Cluster.txt", header=T, sep="|")

# Modify Inputs #
genome_len_db$taxonomy_id <- sub(".*\\|","",sub("_","|",sub("^AC_","AC",sub("^NC_","NC",genome_len_db$taxonomy_id))))
genome_len_db<-arrange(genome_len_db, taxonomy_id)

virusCluster$virusName <- gsub(" ","_",sub(".*\\|","",sub(" ","|",virusCluster$virusName)))

taxon_name <- system("cat */BLAST/*.second_blasat.firstHit.virusOnly.txt | cut -f 12 | sort | uniq | sed -e 's/ /_/g'", intern=T)
genome_len_db <- genome_len_db[genome_len_db$taxonomy_id %in% taxon_name,]

for (virusCluster_num in 1:dim(virusCluster)[1]) {
	genome_len_db$taxonomy_id <- sub(virusCluster[virusCluster_num,1],virusCluster[virusCluster_num,3],genome_len_db$taxonomy_id)
}

genome_len_db <- as.data.table(genome_len_db)
genome_len_db <- as.data.frame(genome_len_db[genome_len_db[, .I[len == max(len)], by=taxonomy_id]$V1])

for ( sample_num in 1:dim(mapped_read)[1] ) {

	sample <- rownames(mapped_read)[sample_num]
	print(sample)

	# Load second Blast Results
	blast <- read.table(paste0(sample,"/BLAST/",sample,".second_blasat.firstHit.virusOnly.txt"), header=F, sep='\t')
	blast[,12] <- gsub(" ", "_",blast[,12])
	dim(blast)
	# Merging Clusters #
	for (virusCluster_num in 1:dim(virusCluster)[1]) {
		blast[,12] <- sub(virusCluster[virusCluster_num,1],virusCluster[virusCluster_num,3],blast[,12])
	}

	# incorrectly paired read remove
	readName <- unique(gsub("/[12]","",blast[,1]))
	
	paired_readName <- c()
	for ( read in readName) {
		if ( paste0(read,"/1") %in% blast[,1] & paste0(read,"/2") %in% blast[,1]) {
			paired_readName <- c(paired_readName, read)  
		}
	}

	blastPairedFiltered <- data.frame()
	
	for (read in readName ){
		if ( read %in% paired_readName ) {
			if ( blast[blast$V1==paste0(read,"/1"),12]==blast[blast$V1==paste0(read,"/2"),12]) { 
				blastPairedFiltered <- rbind(blastPairedFiltered,blast[blast$V1==paste0(read,"/1"),])
				blastPairedFiltered <- rbind(blastPairedFiltered,blast[blast$V1==paste0(read,"/2"),])
			} else {
				print(read)
			}
		} else {
			blastPairedFiltered <- rbind(blastPairedFiltered,blast[blast$V1==paste0(read,"/1"),])
			blastPairedFiltered <- rbind(blastPairedFiltered,blast[blast$V1==paste0(read,"/2"),])
		}
	}
	dim(blastPairedFiltered)
	# Count Filtering  
	blastCount<-as.data.frame(blastPairedFiltered[,12])
	colnames(blastCount) <- "taxonomy_id"
	blastCount <- arrange(blastCount, taxonomy_id)
	blastCount <- count(blastCount)
	blastCount <- blastCount[blastCount$freq > countCutoff,]
	blastPairedCountsFiltered <- blastPairedFiltered[blastPairedFiltered$V12%in%blastCount$taxonomy_id,]
	dim(blastPairedCountsFiltered)

	## binning filtering 	
	blastPairedCountsBinnigFiltered <- data.frame()
	
	blastBedtmp <- blastPairedCountsFiltered[,c(12,6,7)]
	blastBed <- blastBedtmp

	for (i in 1:dim(blastBed)[1]) {
		blastBed[i,1] <- blastBedtmp[i,1]
		blastBed[i,2] <- min(blastBedtmp[i,c(2,3)])
		blastBed[i,3] <- max(blastBedtmp[i,c(2,3)])
	}

	colnames(blastBed) <- c("taxon","start","end")
	blastBed$start <- as.numeric(blastBed$start)
	blastBed$end <- as.numeric(blastBed$end)
	blastBed <- blastBed[c(order(blastBed$taxon, blastBed$start, blastBed$end)),]
	write.table(blastBed, "blast.bed.tmp", quote=FALSE, row.names=F, col.names=F,sep='\t')

	binning <- data.frame()
	for (j in unique(blastBed[,1]) ){
		for (i in 1:10 ) {
			binning <- rbind(binning, c(j, as.integer(genome_len_db[genome_len_db$taxonomy_id==j,2]/10*(i-1)+1), as.integer(genome_len_db[genome_len_db$taxonomy_id==j,2]/10*i)))
		}
	}

	write.table(binning, "binning.tmp", quote=FALSE, row.names=F, col.names=F,sep='\t')

	bedtoolsCoverages_tmp <- system(paste0("bedtools coverage -nonamecheck -a binning.tmp -b blast.bed.tmp"), intern=TRUE)
	bedtoolsCoverages <- data.frame()
	for (i in 1:length(bedtoolsCoverages_tmp) ) {
		bedtoolsCoverages <- rbind(bedtoolsCoverages, unlist(strsplit(bedtoolsCoverages_tmp[i], split='\t')))
	}

	colnames(bedtoolsCoverages) <- c("taxon","start","end","countA","countB","lengthA","fraction")
	binningCount <- count(bedtoolsCoverages[bedtoolsCoverages$countA > 0,]$taxon)
	final_taxonId <- binningCount[binningCount$freq > binningCutoff,]$x
	
	blastPairedCountsBinnigFiltered <- blastPairedCountsFiltered[blastPairedCountsFiltered$V12 %in% final_taxonId,]
	dim(blastPairedCountsBinnigFiltered)

	write.table(blastPairedCountsBinnigFiltered, paste0(sample,"/BLAST/",sample,".second_blast.firstHit.virusOnly.pairedCountBinningFiltered.txt"), col.names=F, sep='\t', row.names=F, quote=F)
	
}
