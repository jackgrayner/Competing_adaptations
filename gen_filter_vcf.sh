#generate vcf from gstacks catalog
populations -P ./ -M original_popmap.tsv -O final_jan23 -r 0.75 -p 2 --min-maf 0.1 --vcf --plink -t 8 

#filter by genotype quality
VCF=populations.snps
vcftools --vcf ./$VCF.vcf --minGQ 20 --recode --recode-INFO-all --out $VCF
VCF1=$VCF.recode

#separate sex/autosomes/header
grep -v '^#' ./$VCF1.vcf | grep -v "^Scaffold_1" > $VCF1.auto.vcf
grep -v '^#' ./$VCF1.vcf | grep "^Scaffold_1" > $VCF1.X.vcf
grep '^#' ./$VCF1.vcf > vcfheader
rm populations.snps.recode.vcfl

#polarise X calls by allelic depth (assuming format is "GT:DP:AD1,AD2:GQ:GL")
Rscript polarise_X.R $VCF1.X.vcf polarised.$VCF1.X.vcf
rm $VCF1.X.vcf

#recombine X/autosomes
cat polarised.$VCF1.X.vcf $VCF1.auto.vcf > polarised.$VCF1.vcf
rm polarised.$VCF1.X.vcf
rm $VCF1.auto.vcf

#filter SNP calls relative to mean sequencing depth (assuming format is "GT:DP:AD1,AD2:GQ:GL")
Rscript filter_meanDepth_sex.R polarised.$VCF1.vcf polarised.$VCF1.depthfilt.vcf1
cat vcfheader polarised.$VCF1.depthfilt.vcf1 | grep "^#\|0/1\|1/1" > polarised.$VCF1.depthfilt.vcf
rm polarised.$VCF1.depthfilt.vcf1

plink --vcf polarised.$VCF1.depthfilt.vcf --make-bed --allow-extra-chr --out $VCF.final 
cp plink_correctedsex.fam ./$VCF.final.fam
plink --bfile $VCF.final --mind 0.6 --geno 0.3 --assoc fisher --allow-extra-chr --adjust --out CwPop_gwas 
