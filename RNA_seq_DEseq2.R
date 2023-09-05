library(DESeq2)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

setwd('/Users/Jack/Desktop/PD/Cw crosses/rna/new_genome/new_anno')
pheno<-read.csv('../phenotype.csv',row.names=1)
pheno$Cw<-factor(substr(pheno$Phenotype,start=6,stop=7))
pheno$Fw<-factor(substr(pheno$Phenotype,start=3,stop=4))

cts<-read.csv('gene_count_matrix.csv',h=T,row.names="gene_id")
#remove genes not expressed at >1 count per mill. in at least 4 samples
cts<-cts[rowSums(cpm(cts)>1)>=4, ]
nrow(cts)

dds <- DESeq(DESeqDataSetFromMatrix(countData = cts,
                              colData = pheno,
                              design= ~ Cw + Fw))

##### QC
vsd<-vst(dds)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("Cw", "Fw"))
####

#### Remove samples S23 & S24 based on aberrant PCA/sample distance correlations
dds <- DESeq(DESeqDataSetFromMatrix(countData = cts[,!colnames(cts) %in% c("S23","S24")],
                              colData = pheno[!rownames(pheno) %in% c("S23","S24"),],
                              design= ~ Cw + Fw))

res.cw <- results(dds, name="Cw_Wt_vs_Cw",alpha=0.05)
resOrdered.cw <- res.cw[order(res.cw$pvalue),]
summary(resOrdered.cw)

res.fw <- results(dds, name="Fw_Nw_vs_Fw",alpha=0.05)
resOrdered.fw <- res.fw[order(res.fw$pvalue),]
summary(resOrdered.fw)


#merge with position info
s.pos<-read.table('stringtie_pos.tsv',h=T)
s.pos<-s.pos[s.pos$gene %in% rownames(cts),]
s.pos<-data.frame(s.pos %>% group_by(chr,gene) %>% summarize(mean = mean(start)))
cw.dge<-data.frame(res.cw)
cw.dge$gene<-rownames(res.cw)
cw.dge.pos<-merge(cw.dge,s.pos,by='gene')
cw.dge.pos<-cw.dge.pos[!is.na(cw.dge.pos$padj),]
cw.dge.pos$DE<-"NS"
cw.dge.pos[cw.dge.pos$padj<0.05,]$DE<-"S1"
cw.dge.pos$scaffold<-factor(sub("Scaffold_","",cw.dge.pos$chr),
                            levels=c(1:14))

fw.dge<-data.frame(res.fw)
fw.dge$gene<-rownames(res.fw)
fw.dge.pos<-merge(fw.dge,s.pos,by='gene')
fw.dge.pos<-fw.dge.pos[!is.na(fw.dge.pos$padj),]
fw.dge.pos$DE<-"NS"
fw.dge.pos[fw.dge.pos$padj<0.05,]$DE<-"S1"
fw.dge.pos$scaffold<-factor(sub("Scaffold_","",fw.dge.pos$chr),
                            levels=c(1:14))

#plot DE genes across first 14 scaffolds
levels(cw.dge.pos$scaffold)[1]<-"1 (X)"
levels(fw.dge.pos$scaffold)[1]<-"1 (X)"

cw.twas_hor<-ggplot(cw.dge.pos,aes(x=mean/1e+06,y=-log10(padj),colour=DE))+
  geom_point(alpha=1,size=0.75)+theme_bw()+theme(panel.grid=element_blank(),
                                                   axis.ticks=element_blank(),
                                                   axis.text.x=element_blank(),
                                                   legend.position='none')+
  facet_grid(.~scaffold,scales='free_x',space='free_x')+
  xlab('Genome position')+ylab('-log10 (Padj.)')+
  scale_colour_manual(values=c('#555555','red','#603c75'))+
  geom_rug(data=cw.dge.pos[!cw.dge.pos$DE=="NS",],aes(x=mean/1e+06),
           sides='b',alpha=0.3)

fw.twas_hor<-ggplot(fw.dge.pos,aes(x=mean/1e+06,y=-log10(padj)))+
  geom_point(alpha=1,size=0.75,aes(colour=DE))+theme_bw()+
  theme(panel.grid=element_blank(),axis.text.x=element_blank(),
                                                   axis.ticks=element_blank(),
                                                   legend.position='none')+
  facet_grid(.~scaffold,scales='free_x',space='free_x')+
  xlab('Genome position')+ylab('-log10 (Padj.)')+
  scale_colour_manual(values=c('#555555','red'))+
  geom_rug(data=fw.dge.pos[!fw.dge.pos$DE=="NS",],aes(x=mean/1e+06),
           sides='b',alpha=0.3,colour='red')

grid.arrange(cw.twas_hor+ggtitle('Curly-wing DEGs'),
             fw.twas_hor+ggtitle('Flatwing DEGs'),nrow=2)

ggsave('CwRAD_DEGs_pos_hor.png',dpi=600,height=9,width=9,
    grid.arrange(cw_all+ggtitle('Curly-wing genetic association'),
                 cw.twas_hor+ggtitle('Curly-wing DEGs'),
                 fw.twas_hor+ggtitle('Flatwing DEGs'),nrow=3)
)


#DEcw heatmap
ntd <- normTransform(dds)
select<-which(res.cw$padj<0.05)
df <- as.data.frame(colData(dds)[,c("Cw","Fw")])
pheatmap(assay(ntd)[select,],cluster_rows=TRUE,show_rownames=FALSE,
         cluster_cols=TRUE,annotation_col=df,scale='row')

#DEfw heatmap
ntd <- normTransform(dds)
select<-which(res.fw$padj<0.05)
df <- as.data.frame(colData(dds)[,c("Cw","Fw")])
pheatmap(assay(ntd)[select,],cluster_rows=TRUE,show_rownames=FALSE,
         cluster_cols=TRUE,annotation_col=df,scale='row')
