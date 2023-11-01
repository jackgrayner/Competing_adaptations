vcf1 = commandArgs(trailingOnly=TRUE)[1]
vcf2 = commandArgs(trailingOnly=TRUE)[2]

vcf<-read.table(vcf1,h=F)
#find AD field
#format "GT:DP:AD1,AD2:GQ:GL"

#keep only biallelic
vcf<-vcf[nchar(vcf[,5])==1,]

#identify male samples
sexes<-read.table('plink_correctedsex.fam',h=F)
males<-which(sexes[,5]=="1")

#loop to remove genotypes for samples with heterozygote genotypes
for (SNP in c(1:nrow(vcf))){
  #the first 9 columns are not samples
  #read genotype from respective locus,sample
  for (sample in (males+9)){
    if(grepl("0/1",vcf[SNP,sample])){
      gen=vcf[SNP,sample]
      #extract genotype and allelic imbalance P-value
      #heterozygous
      ind1 = unlist(gregexpr(pattern = ":", text = gen))
      AD1<-substr(gen,ind1[length(ind1)-2] + 1, ind1[length(ind1)-1]-1)
      Fref<-as.numeric(substr(AD1, 0,unlist(gregexpr(pattern=",",AD1))-1))
      Falt<-as.numeric(substr(AD1, unlist(gregexpr(pattern=",",AD1))+1, nchar(AD1)))
      if (binom.test(c(Fref,Falt))$p.value>0.05) {
        #if no sig. allele imbalance, genotype is null
        vcf[SNP,sample]<-sub("0/1", "./.", vcf[SNP,sample])
        #if sig. allele imbalance, rewrite genotype as homozygous based on major allele
      } else {
        #extract specific allelic depths
        if (Fref>Falt){
          vcf[SNP,sample]<-sub("0/1", "0/0", vcf[SNP,sample])
        } else {
          vcf[SNP,sample]<-sub("0/1", "1/1", vcf[SNP,sample])
        }
      }
    }
  }
}

write.table(vcf,file=vcf2,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
