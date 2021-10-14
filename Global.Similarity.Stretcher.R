# Env
library(data.table)

# Args
args = commandArgs(trailingOnly=TRUE)
Global_start <- as.numeric(args[1])
Global_end <- as.numeric(args[2])
virus_fna <- args[3]

# Data Load
Allref <- fread(virus_fna, header=F)
HeaderIndex <- grep("^>",Allref$V1)
HeaderIndex <- c(HeaderIndex, dim(Allref)[1]+1)
HeaderIndex_Node <- HeaderIndex[Global_start:Global_end]

similarity_db <- data.frame()

# calculate global similarity
for ( i in 1:as.numeric(length(HeaderIndex_Node)-1)) {
  seq1_start <- HeaderIndex_Node[i]
  seq1_end <- HeaderIndex_Node[i+1]-1
  seq1 <- Allref[c(seq1_start:seq1_end),]
  write.table(as.data.frame(seq1), paste0("seq1_",Global_start,"_",Global_end,".tmp"), sep='\t', quote=F, col.names=F, row.names=F)
  for ( j in match(HeaderIndex_Node[i], HeaderIndex):as.numeric(length(HeaderIndex)-1)) {
    print(paste0(Allref[HeaderIndex_Node[i]], Allref[HeaderIndex[j]]))
    seq2_start <- HeaderIndex[j]
    seq2_end <- HeaderIndex[j+1]-1
    seq2 <- Allref[c(seq2_start:seq2_end),]
    write.table(as.data.frame(seq2),paste0("seq2_",Global_start,"_",Global_end,".tmp") , sep='\t', quote=F, col.names=F, row.names=F)
    
    system(paste0("stretcher seq1_",Global_start,"_",Global_end,".tmp"," seq2_",Global_start,"_",Global_end,".tmp -outfile stretcher.result.",Global_start,"_",Global_end,".tmp"))
    similarity <- system(paste0("grep Similarity ","stretcher.result.",Global_start,"_",Global_end,".tmp"," | sed -e 's/.*(//g' | sed -e 's/%)//g'"), intern=T)
    
    
    write(paste(as.character(Allref[HeaderIndex[i],]), as.character(Allref[HeaderIndex[j],]), as.character(similarity), sep="|"), file=paste0("virus.global.similarity.virusNum",Global_start,"_",Global_end,".txt") ,append=TRUE)
  }
}