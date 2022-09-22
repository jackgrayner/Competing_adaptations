vcf<-read.table('pop.snps.prefilter.vcf',h=T)

#find depth field
extr.depth<-function(strings){
  ind = unlist(gregexpr(pattern = ":", text = strings))
  if (length(ind) < 4){NA}
  else{substr(strings, ind[length(ind) - 3] + 1, ind[length(ind)-2] - 1)}
}
means<-read.table('averageList.txt',h=F)[,1]
#means<-rnorm(380,mean=30,sd=5)

sexes<-read.table('plink.fam',h=F)[,5]

#loop to remove genotypes for samples with abnormal sequencing depth
for (rowNum in c(1:nrow(vcf))){
  #the first 9 columns are not samples
  for (sampNum in c(10:ncol(vcf))) {
    #-9 to get to sample ID
    d=sampNum-9
    #read genotype from respective locus,sample
    gen=vcf[rowNum,sampNum]
    #extract mean depth of respective sample
    meanDepth<-means[d]
    #if locus is NOT null
    if (!is.na(as.numeric(sapply(gen,extr.depth)))){
      #if locus is not X-linked, or if sex is female
      if (!vcf[rowNum,1]=="Scaffold_1" | sexes[d]=="2") {
        #if depth at this locus is > 3x of sample mean, replace as null
        if (as.numeric(sapply(gen,extr.depth)) > (3 * meanDepth)) {
          vcf[rowNum,sampNum]<-"./."
        }
        #else if depth at this locus is < 1/3x of sample mean, replace as null
        else if (as.numeric(sapply(gen,extr.depth)) < (1/3 * meanDepth)) {
          vcf[rowNum,sampNum]<-"./."
        }
        #else if depth <10, replace as null
        else if (as.numeric(sapply(gen,extr.depth)) < 10) {
          vcf[rowNum,sampNum]<-"./."
        }
        #else if locus is X-linked and sample is male, halve filtering criteria
      } else if (vcf[rowNum,1]=="Scaffold_1" & sexes[d]=="1") {
        #if X-linked, alter criteria
        #if depth at this locus is > 1.5x of sample mean, replace as null
        if (as.numeric(sapply(gen,extr.depth)) > (1.5 * meanDepth)) {
          vcf[rowNum,sampNum]<-"./."
        }
        #else if depth at this locus is < 1/6x of sample mean, replace as null
        else if (as.numeric(sapply(gen,extr.depth)) < (1/6 * meanDepth)) {
          vcf[rowNum,sampNum]<-"./."
        }
        #else if depth <10, replace as null
        else if (as.numeric(sapply(gen,extr.depth)) < 10) {
          vcf[rowNum,sampNum]<-"./."
        }
      }
    }
  }
  if (rowNum %% 1000==0){
    print(paste(round(rowNum/(nrow(vcf))*100,digits=2),"% done"))
  }
}


write.table(vcf,file='pop.snps.filtered.vcf',sep='\t',row.names=FALSE,quote=FALSE)
