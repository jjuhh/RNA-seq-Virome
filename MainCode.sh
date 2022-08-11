# Human Virome analysis pipeline #
## Kumata, Ito and Sato in prep ##
#
# input : RNA-seq data
# output : viral read count

##requirements##
# Require the following programs
#   trimmomatic
#   STAR
#   featureCount
#   samtools
#   picard
#   python2
#   blastn
# Requiew the following genome
#   Human genome
#   prokaryota genome

# set directory #
sample=$1
wd="."
path_to_input=${wd}/${sample}"/unmapped"
path_to_STAR=${wd}/${sample}"/Align"
path_to_unmapped=${wd}/${sample}"/unmapped"
path_to_database="/gmi-l1/_90.User_Data/juhyunk/Files/BlastDB/virus_only_ref/virus.virome.DB"
path_to_viral_hit_fasta=${wd}/${sample}"/BLAST"
path_to_res_first_blast=${wd}/${sample}"/BLAST"
path_to_res_second_blast=${wd}/${sample}"/BLAST"
path_to_count=${wd}/${sample}"/BLAST"
script="/gmi-l1/_90.User_Data/juhyunk/codes/virome/sato/TheSatoLab-Human_Virome_analysis-5140486/code"

# set threshold
threshold=1.0e-10
mkdir ${path_to_res_first_blast}

# 4. First blast search #
blastn -db ${path_to_database}/virus_genome_for_first_blast.db -query ${path_to_unmapped}/${sample}.unmapped.fasta -out ${path_to_res_first_blast}/${sample}.first_blast.txt -word_size 11 -outfmt 6 -num_threads 50 -evalue ${threshold}

# 5. read. extract
cut -f 1 ${path_to_res_first_blast}/${sample}.first_blast.txt > ${path_to_res_first_blast}/firstBlast.readName.txt
grep -A 1 -f ${path_to_res_first_blast}/firstBlast.readName.txt  ${path_to_unmapped}/${sample}.unmapped.fasta | grep -G -v "^-" > ${path_to_res_first_blast}/${sample}.first_blast_hit.fasta
#rm ${path_to_res_first_blast}/firstBlast.readName.txt

# 6. second blast search
blastn -db ${path_to_database}/virus_genome_for_second_blast.db -query ${path_to_res_first_blast}/${sample}.first_blast_hit.fasta -out ${path_to_res_first_blast}/${sample}.second_blast.txt -word_size 11 -outfmt '6 qseqid sseqid pident evalue bitscore sstart send staxids sscinames scomnames sskingdoms stitle' -num_threads 50 -evalue ${threshold}
