rm(list=ls())
####PCA
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
