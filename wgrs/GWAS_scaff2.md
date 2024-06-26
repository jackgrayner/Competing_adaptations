```sh
vcf=$1
ref=~/scratch/Cricket.curated.scaff_v2.fasta
bamlist=./bamlist.txt
bcftools mpileup -Ov -f $ref -d 100000 --annotate FORMAT/AD,FORMAT/DP -b $bamlist -r Scaffold_2 | bcftools call -m -v -Oz -f GQ -o $vcf.vcf.gz
vcftools --gzvcf $vcf.vcf.gz --maf 0.1 --minDP 5 --maxDP 100 --minQ 20 --minGQ 20 --recode --recode-INFO-all --out $vcf
bgzip $vcf.recode.vcf
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' $vcf.recode.vcf.gz -O z -o $vcf.recode_anno.vcf.gz
bcftools index $vcf.recode_anno.vcf.gz
```
```sh
plink --allow-extra-chr --vcf $vcf.recode_anno.vcf.gz \
--chr Scaffold_2 --double-id --thin 0.001 --make-bed --snps-only --allow-no-sex \
--out $vcf.thin --export HV
java -jar ~/scratch/apps/Haploview.jar -memory 5000
```

```sh
bcftools merge -r Scaffold_2 -m id ../popgenHUH_scaff2.recode_anno.vcf.gz ../popgenKCG_scaff2.recode_anno.vcf.gz ../popgenOCC_scaff2.recode_anno.vcf.gz > HUH_KCG_OCC.vcf.gz
cat ../popgenHUH_scaff2.fam ../popgenKCG_scaff2.fam ../popgenOCC_scaff2.fam > HUH_KCG_OCC_merged.fam
```

```sh
plink --vcf HUH_KCG_OCC.vcf.gz --geno 0.1 ---chr Scaffold_2 --maf 0.15 --double-id --make-bed --allow-extra-chr --out HUH_KCG_OCC
cp ./cw.fam HUH_KCG_OCC.fam
gemma -lm 2 -c Cw_pop_covar.txt -miss 0.1 -bfile HUH_KCG_OCC -o HUH_KCG_OCC
#this uses likelihood ratio test LM with a covariate of population
```

```sh
plink --vcf HUH_KCG_OCC.vcf.gz --double-id --allow-extra-chr \
--geno 0.2 --maf 0.15 --chr Scaffold_2 --from-mb 7.5 --to-mb 80 \
--make-bed --pca --out HUH_KCG_OCC
```

```sh
#plot heterozygosity
populations -V ./$vcf.recode_anno.vcf.gz -O ./ -R 0.5 -M $popmap
```