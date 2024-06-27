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

```R
library(ggplot2)
dat<-read.table('HUH_KCG_OCC.assoc.txt',h=T)
dat<-dat[nchar(dat$allele1)==1 & nchar(dat$allele0)==1,]
dat<-dat[dat$af>0.15,]

dat$Padj<-p.adjust(dat$p_lrt,method="BH")
dat$outlier<-"N"
dat[dat$Padj<quantile(dat$Padj,0.001),]$outlier<-"Y"

serps<-read.table('serpins_GO0004867.tsv',h=T)
serps.2<-serps[serps$scaffold=="Scaffold_2",]

allpops.gwas<-ggplot(dat,aes(x=ps/1e+06,y=-log10(p_lrt),colour=sig))+
	geom_rect(xmin=7.5,xmax=80,ymin=(-1),ymax=10,fill='#e9f5dc',alpha=0.3,colour='white')+
	geom_rug(data=serps.2,aes(x=start/1e+06),inherit.aes=FALSE,colour='black',sides='top')+
	geom_point(size=0.5,alpha=0.5)+theme_bw()+
	scale_colour_manual(values=c('#666666','#b0473f'))+
	theme(legend.position='none',panel.grid=element_blank())+
	scale_x_continuous(breaks=seq(0,250,by=20))+
	xlab('Chr2 pos. (MB)')+ylab('-log10(P)')

```

```sh
#plot heterozygosity
populations -V ./$vcf.recode_anno.vcf.gz -O ./ -R 0.5 -M $popmap
```
