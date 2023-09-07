#!/bin/bash
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

FILES="./*.fq.gz"
for f in $FILES
do
 bwa mem -M -t 16 ~/scratch/Cricket.curated.scaff_v2.fasta $f > $f.sam
done
