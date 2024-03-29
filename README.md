# RNA-seq-Virome

## Introduction

A viral genomic sequence database was constructed as follows.    
1. Of the viral sequences, sequence regions that highly resemble that human or bacterial genomes were masked (i.e., replaced by the sequence "NNN....").   
2. redundant viral sequence (i.e., sequences disclosing >90% global sequence identity with each other) were concatenated.    

1st Blastn DB : viral genomes prepared above.       
2nd Blastn DB : viral genomes prepared above + human genome(CHM13) + bacterial genome.           
ps. 2021 10 13.          
CHM13 v.1.1 used, but that has not chrY.          

## Masking Virus reference fasta file
1. For masking, to find resemble region between bacterial or human and viral genome, conduct BLASTn.       
The sequence regions to be masked were determined by a local sequence similarity search using BLASTn (ver 2.11.0+)            
Parameters : word size(11), E value(1.0e-03)               
Sources : human genome(CHM13 v.1.1), bacterial genome(total ), viral genome(total )#            
2. Formatting to BED and merge the resemble regions.         
Use bedtools merge command.         
Before use bedtools, All of contig and bp column should be sorted.              
Input : Blastn outfmt 6 files           
Output : bed file to mask          
3. masking fasta file         
Input : fasta file, bed file prepared above            
Output : fasta file masked     

### Find resemble region between virus and human CHM13 reference 
```
blastn -query virus_genome.fna \ 
-db chm13.draft_v1.1.fasta  \
-out Human.virus.resemble.regions.identification \
-word_size 11 -outfmt '6 qseqid sseqid pident evalue bitscore qstart qend sstart send staxids sscinames scomnames sskingdoms stitle' \
-num_threads 30 -evalue 1e-3 > Human.virus.resemble.regions.identification.stout 2> Human.virus.resemble.regions.identification.sterr
```
###  Formatting to BED and merge the resemble regions.   
```
bedtools merge -i <(cut -f 1,6,7 Human.virus.resemble.regions.identification| sort -k1,1V -k2n ) > Human.virus.resemble.regions.identification.merge.bed
```
### Masking resemble region
```
bedtools maskfasta \
-fi virus_genome.fna \
-bed Human.virus.resemble.regions.identification.merge.bed \
-fo virus_genome.masking.fna \
-fullHeader
```
### Find resemble region between virus and human bacterial reference 
```
blastn \
-query virus_genome.fna \
-db ref_prok_rep_genomes \
-out Human.virus.prokaryotes.resemble.regions.identification \
-word_size 11 \
-outfmt '6 qseqid sseqid pident evalue bitscore qstart qend sstart send staxids sscinames scomnames sskingdoms stitle' \
-num_threads 30 \
-evalue 1e-3 > Human.virus.prokaryotes.resemble.regions.identification.stout 2> Human.virus.prokarytoes.resemble.regions.identification.sterr
```
###  Formatting to BED and merge the resemble regions.   
```
bedtools merge -i <(cut -f 1,6,7 Human.virus.prokaryotes.resemble.regions.identification| sort -k1,1V -k2n ) > Human.virus.prokaryotes.resemble.regions.identification.merge.bed
```
### Masking resemble region
```
bedtools maskfasta \
-fi virus_genome.making.fna \
-bed Human.virus.prokaryotes.resemble.regions.identification.merge.bed \
-fo virus_genome.Human_prokaryotes_masking.fna \
-fullHeader
```
## Blast search
``` 
mkdir ${path_to_res_first_blast}
blastn -db ${path_to_database}/virus_genome_for_first_blast.db -query ${path_to_unmapped}/${sample}.unmapped.fasta -out ${path_to_res_first_blast}/${sample}.first_blast.txt -word_size 11 -outfmt 6 -num_threads 50 -evalue ${threshold}
```
## Similarity Calculation between virus detected by BLAST 
```
In debugging
```

## First Hit read Extreact
```
fgrep -A 1 -f <(cut -f 1 ${path_to_res_first_blast}/${sample}.first_blast.txt) ${path_to_unmapped}/${sample}.unmapped.fasta | grep -G -v "^-" > ${path_to_res_first_blast}/${sample}.first_blast_hit.fasta
```
# 6. second blast search
```
blastn -db ${path_to_database}/virus_genome_for_second_blast.db -query ${path_to_res_first_blast}/${sample}.first_blast_hit.fasta -out ${path_to_res_first_blast}/${sample}.second_blast.txt -word_size 11 -outfmt '6 qseqid sseqid pident evalue bitscore sstart send staxids sscinames scomnames sskingdoms stitle' -num_threads 50 -evalue ${threshold}
```
