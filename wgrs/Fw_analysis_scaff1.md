### Identify, call and filter variants
```sh
vcf=$1#provide name of vcf
ref=~/scratch/Cricket.curated.scaff_v2.fasta
bamlist=./bamlist.txt#list of bam files
bcftools mpileup -Ov -f $ref -d 100000 --annotate FORMAT/AD,FORMAT/DP -b $bamlist -r Scaffold_1 | bcftools call --ploidy 1 -m -v -Oz -f GQ -o $vcf.vcf.gz
vcftools --gzvcf $vcf.vcf.gz --maf 0.1 --minDP 5 --maxDP 70 --minQ 20 --minGQ 20 --recode --recode-INFO-all --out $vcf
bgzip $vcf.recode.vcf
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' $vcf.recode.vcf.gz -O z -o $vcf.recode_anno.vcf.gz
```

###  View LD
```sh
plink --allow-extra-chr --vcf $vcf.recode_anno.vcf.gz \
--chr Scaffold_1 --double-id --remove remove.txt --thin 0.001 --make-bed --snps-only --geno 0.5 --allow-no-sex \
--out $vcf.thin --export HV
java -jar ~/scratch/apps/Haploview.jar -memory 5000
```

### Association test
```sh
cp ./fw.fam $vcf.fam
gemma -lm 2 -miss 0.15 -bfile $vcf -o $vcf
```
### Plot association test results
```R
kcg.fw<-read.table('popgenKCG_X.assoc.assoc.txt',h=T)
kcg.fw<-kcg.fw[nchar(kcg.fw$allele1)==1 & nchar(kcg.fw$allele0)==1,]
kcg.fw<-kcg.fw[kcg.fw$af>0.2,]
kcg.fw$Padj<-p.adjust(kcg.fw$p_lrt,method="BH")
kcg.fw$outlier<-"N"
kcg.fw[kcg.fw$Padj<quantile(kcg.fw$Padj,0.001),]$outlier<-"Y"

#add extra DPs for KCG plot
scaleFUN <- function(x) sprintf("%.2f", x)

gwas.fw.kcg<-ggplot(kcg.fw,aes(x=ps/1e+06,y=-log10(p_lrt),colour=sig))+
	geom_rect(xmin=95.7,xmax=253.2 ,ymin=(-1),ymax=12,fill='#e9f5dc',alpha=0.3,colour='white')+
	geom_vline(xintercept=259.2,colour='orange',size=1,alpha=1.5)+
	geom_point(size=0.7,alpha=0.5)+theme_bw()+
	scale_colour_manual(values=c('#666666','#b0473f'))+
	theme(legend.position='none',panel.grid=element_blank())+
	scale_x_continuous(breaks=seq(0,350,by=20))+
	xlab('ChrX pos. (MB)')+ylab('-log10(P)')+scale_y_continuous(labels=scaleFUN)

occ.fw<-read.table('popgenOCC_X.assoc.assoc.txt',h=T)
occ.fw<-occ.fw[nchar(occ.fw$allele1)==1 & nchar(occ.fw$allele0)==1,]
occ.fw<-occ.fw[occ.fw$af>0.2,]
occ.fw$Padj<-p.adjust(occ.fw$p_lrt,method="BH")
occ.fw$outlier<-"N"
occ.fw[occ.fw$Padj<quantile(occ.fw$Padj,0.001),]$outlier<-"Y"

gwas.fw.occ<-ggplot(occ.fw,aes(x=ps/1e+06,y=-log10(p_lrt),colour=sig))+
	geom_rect(xmin=95.7,xmax=253.2 ,ymin=(-1),ymax=12,fill='#e9f5dc',alpha=0.3,colour='white')+
	geom_vline(xintercept=259.2,colour='orange',size=1,alpha=1.5)+
	geom_point(size=0.7,alpha=0.5)+theme_bw()+
	scale_colour_manual(values=c('#666666','#b0473f'))+
	theme(legend.position='none',panel.grid=element_blank())+
	scale_x_continuous(breaks=seq(0,350,by=20))+
	xlab('ChrX pos. (MB)')+ylab('-log10(P)')

```
