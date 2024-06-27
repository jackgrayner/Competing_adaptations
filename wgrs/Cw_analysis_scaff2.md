### Identify, call and filter variants
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

### View linkage
```sh
plink --allow-extra-chr --vcf $vcf.recode_anno.vcf.gz \
--chr Scaffold_2 --double-id --thin 0.001 --make-bed --snps-only --allow-no-sex \
--out $vcf.thin --export HV
java -jar ~/scratch/apps/Haploview.jar -memory 5000
```

### Merge variant files
```sh
bcftools merge -r Scaffold_2 -m id ../popgenHUH_scaff2.recode_anno.vcf.gz ../popgenKCG_scaff2.recode_anno.vcf.gz ../popgenOCC_scaff2.recode_anno.vcf.gz > HUH_KCG_OCC.vcf.gz
cat ../popgenHUH_scaff2.fam ../popgenKCG_scaff2.fam ../popgenOCC_scaff2.fam > HUH_KCG_OCC_merged.fam
```

### Run association tests
```sh
plink --vcf HUH_KCG_OCC.vcf.gz --geno 0.1 ---chr Scaffold_2 --maf 0.15 --double-id --make-bed --allow-extra-chr --out HUH_KCG_OCC
cp ./cw.fam HUH_KCG_OCC.fam
gemma -lm 2 -c Cw_pop_covar.txt -miss 0.1 -bfile HUH_KCG_OCC -o HUH_KCG_OCC
#this uses likelihood ratio test LM with a covariate of population
```

### Run PCA
```sh
plink --vcf HUH_KCG_OCC.vcf.gz --double-id --allow-extra-chr \
--geno 0.2 --maf 0.15 --chr Scaffold_2 --from-mb 7.5 --to-mb 80 \
--make-bed --pca --out HUH_KCG_OCC
```

### Plot association test results
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

### Plot PCA
```R
library(tidyr)
library(tidyverse)
library(ggbeeswarm)
library(ggridges)

pca <- read_table2("HUH_KCG_OCC.eigenvec", col_names = FALSE)
eigenval <- scan("/HUH_KCG_OCC.eigenval")
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

labs<-read.table("/merged_HUH_KCG_OCC_fam.txt")
pca$labs<-labs[,6]
pca$labs<-factor(pca$labs)
levels(pca$labs)<-c("Wt","Cw")
pca$pop<-c(rep("Hawaii.UH",30),rep("Kauai.CG",30),rep("Oahu.CC",30))

pops.pc1<-ggplot(pca, aes(x=PC1, y=pop,fill=labs))+
  scale_fill_manual(values=c('#C06C84','#355C7D'))+
  geom_density_ridges(alpha=0.5,scale=0.7,
                      point_size = 1, point_alpha = 1,aes(point_colour=labs),
                      jittered_points = TRUE,colour='white',size=0)+theme_minimal()+
  theme(panel.grid=element_blank(),axis.title.y=element_blank(),
        legend.title=element_blank(),legend.position='left')+
  scale_colour_manual(values=c('#C06C84','#355C7D'))+
  theme(plot.margin=margin(c(5.5,50,5.5,0),'points'))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  scale_discrete_manual("point_colour",values=c('#C06C84','#355C7D'))

ggsave('merged_PC1.png',dpi=300,height=5,width=4,
       plot=pops.pc1)
```
### plot heterozygosity
```sh
populations -V ./$vcf.recode_anno.vcf.gz -O ./ -R 0.5 -M $popmap
```
