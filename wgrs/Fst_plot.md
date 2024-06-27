```sh
vcftools --gzvcf HUH_KCG_OCC.vcf.gz
	--chr Scaffold_2
	--fst-window-size 10000
	--fst-window-step 10000
	--weir-fst-pop HUH_inv.txt
	--weir-fst-pop KCG_inv.txt
	--keep HUH_inv.txt
	--keep KCG_inv.txt
	--out HUHinv_KCG_inv.FST
```
```R
fst<-read.table("./HUHinv_KCG_inv.FST.windowed.weir.fst",h=T)
serps<-read.table("./serpins.txt",h=T)
serps.2<-serps[serps$scaffold=="Scaffold_2",]

fst$outlier<-"N"
fst[fst$WEIGHTED_FST>quantile(fst$WEIGHTED_FST,(99.9/100)),]$outlier<-"Y"
fst[fst$WEIGHTED_FST<0,]$WEIGHTED_FST<-0
g.fst<-ggplot(fst,aes(x=BIN_START/1e+06,y=WEIGHTED_FST))+
  geom_rect(xmin=7.5,xmax=80,ymin=(-0.1),ymax=1.2,fill='#def2c7',alpha=0.3,colour='white')+
  geom_rug(data=serps.2[serps.2$start<100e+06,],aes(x=start/1e+06),inherit.aes=FALSE,colour='black',sides='top')+
  geom_point(alpha=0.5,size=1,colour='#555555')+theme_bw()+theme(panel.grid=element_blank())+
  geom_smooth(method="loess",span=0.1,se=TRUE,colour='darkred',size=0.75)+
  ylim(c(0,1))+scale_x_continuous(breaks=seq(0,100,10),limits=c(0,100))+
  xlab("Window start (Mb)")+ylab("Weighted FST")

ggsave('./HUH_KCG.fst.png',plot=g.fst,dpi=600,height=2.75,width=6)
```
