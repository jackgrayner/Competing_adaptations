source activate bwa-mem2
bwa-mem2 mem ~/scratch/Cricket.curated.scaff_v2.fasta $1_merged_R1.fastq.gz $1_merged_R2.fastq.gz | samtools sort -o ./$1.bam

samtools index $1.bam

samtools view -b $1 Scaffold_2 > Chr2_$1.bam
java -jar ~/scratch/apps/picard.jar MarkDuplicates \
      I=Chr2_$1.bam \
      O=MD_Chr2_$1.bam \
      M=MD_Chr2_$1.txt

samtools view -b $1 Scaffold_1 > ChrX_$1.bam
java -jar ~/scratch/apps/picard.jar MarkDuplicates \
      I=ChrX_$1.bam \
      O=MD_ChrX_$1.bam \
      M=MD_ChrX_$1.txt
