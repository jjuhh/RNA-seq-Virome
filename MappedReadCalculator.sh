
for sample in `ls | grep ^SRR`
do echo $sample
mapped_read="`samtools flagstat -@ 20 $sample/Align/$sample*.sorted.dp.bam | sed -n '5p' | cut -d ' ' -f 1`"
echo -e "$sample\t$mapped_read" 
echo -e "$sample\t$mapped_read" >> All.sample.mapped_read
done
