#!/bin/bash
#SBATCH --job-name=bcftools_call
#SBATCH --export=ALL
#SBATCH --mem=16G
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --output=bcftools_call.out
vcf=$1
ref=~/scratch/Cricket.curated.scaff_v2.fasta
bamlist=./bamlist.txt
bcftools mpileup -Ov -f $ref -d 100000 --annotate FORMAT/AD,FORMAT/DP -b $bamlist -r Scaffold_2 | bcftools call -m -v -Ov -f GQ -o $vcf.vcf
vcftools --gzvcf $vcf.vcf.gz --maf 0.1 --minDP 5 --maxDP 100 --minQ 20 --minGQ 20 --recode --recode-INFO-all --out $vcf
bgzip $vcf.recode.vcf
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' $vcf.recode.vcf.gz -O z -o $vcf.recode_anno.vcf.gz

plink --vcf $vcf.recode_anno.vcf.gz --double-id --make-bed --allow-extra-chr --out $vcf
cp cw.fam ./$vcf.fam
plink --bfile $vcf --assoc fisher --double-id --allow-extra-chr --maf 0.15 --geno 0.2 --out $vcf

#LD
plink --allow-extra-chr --vcf $vcf.recode_anno.vcf.gz \
--chr Scaffold_2 --double-id --remove remove.txt --thin 0.001 --make-bed --snps-only --geno 0.5 --allow-no-sex \
--out $vcf.thin --export HV
java -jar ~/scratch/apps/Haploview.jar -memory 5000
