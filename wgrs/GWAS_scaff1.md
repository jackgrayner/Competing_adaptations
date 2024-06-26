```sh
vcf=$1#provide name of vcf
ref=~/scratch/Cricket.curated.scaff_v2.fasta
bamlist=./bamlist.txt#list of bam files
bcftools mpileup -Ov -f $ref -d 100000 --annotate FORMAT/AD,FORMAT/DP -b $bamlist -r Scaffold_1 | bcftools call --ploidy 1 -m -v -Oz -f GQ -o $vcf.vcf.gz
vcftools --gzvcf $vcf.vcf.gz --maf 0.1 --minDP 5 --maxDP 70 --minQ 20 --minGQ 20 --recode --recode-INFO-all --out $vcf
bgzip $vcf.recode.vcf
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' $vcf.recode.vcf.gz -O z -o $vcf.recode_anno.vcf.gz
```

```sh
plink --allow-extra-chr --vcf $vcf.recode_anno.vcf.gz \
--chr Scaffold_1 --double-id --remove remove.txt --thin 0.001 --make-bed --snps-only --geno 0.5 --allow-no-sex \
--out $vcf.thin --export HV
java -jar ~/scratch/apps/Haploview.jar -memory 5000
```

```sh
cp ./fw.fam $vcf.fam
gemma -lm 2 -miss 0.15 -bfile $vcf -o $vcf
```