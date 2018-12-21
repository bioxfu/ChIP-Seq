library(Gviz)
library(GenomicFeatures)
library(RColorBrewer)
col_set <- brewer.pal(9, 'Set1')
options(ucscChromosomeNames=FALSE)

txdb <- loadDb('/home/xfu/Gmatic6/genome/PhytozomeV12/Bdistachyon/annotation/Bdi_txdb.sqlite')
grtrack <- GeneRegionTrack(txdb, collapseTranscripts='longest', transcriptAnnotation='gene', name='Gene', showTitle=F, col='black', fill='black')
gtrack <- GenomeAxisTrack(scale=0.1, col='black')

dtrack1.1 <- DataTrack(range='/home/xfu/Project/CZQ/MeDIP/track/A.bedgraph', name='A', col.mountain='white', fill.mountain=c('white', col_set[1]))
dtrack1.2 <- DataTrack(range='/home/xfu/Project/CZQ/MeDIP/track/C.bedgraph', name='C', col.mountain='white', fill.mountain=c('white', col_set[1]))
dtrack1.3 <- DataTrack(range='/home/xfu/Project/CZQ/MeDIP/track/E.bedgraph', name='E', col.mountain='white', fill.mountain=c('white', col_set[1]))
dtrack1.4 <- DataTrack(range='/home/xfu/Project/CZQ/MeDIP/track/D.bedgraph', name='D', col.mountain='white', fill.mountain=c('white', col_set[1]))

dtrack2.1 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/A_H3K4me3.bedgraph', name='A', col.mountain='white', fill.mountain=c('white', col_set[2]))
dtrack2.2 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/C_H3K4me3.bedgraph', name='C', col.mountain='white', fill.mountain=c('white', col_set[2]))
dtrack2.3 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/E_H3K4me3.bedgraph', name='E', col.mountain='white', fill.mountain=c('white', col_set[2]))
dtrack2.4 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/D_H3K4me3.bedgraph', name='D', col.mountain='white', fill.mountain=c('white', col_set[2]))

dtrack3.1 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/A_H3K9me3.bedgraph', name='A', col.mountain='white', fill.mountain=c('white', col_set[3]))
dtrack3.2 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/C_H3K9me3.bedgraph', name='C', col.mountain='white', fill.mountain=c('white', col_set[3]))
dtrack3.3 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/E_H3K9me3.bedgraph', name='E', col.mountain='white', fill.mountain=c('white', col_set[3]))
dtrack3.4 <- DataTrack(range='/home/xfu/Project/CZQ/ChIP-Seq/track/D_H3K9me3.bedgraph', name='D', col.mountain='white', fill.mountain=c('white', col_set[3]))

dtrack4.1 <- DataTrack(range='/home/xfu/Project/CZQ/RNA-Seq/track/A2-1.bedgraph', name='A', col.mountain='white', fill.mountain=c('white', col_set[4]))
dtrack4.2 <- DataTrack(range='/home/xfu/Project/CZQ/RNA-Seq/track/C1-1.bedgraph', name='C', col.mountain='white', fill.mountain=c('white', col_set[4]))
dtrack4.3 <- DataTrack(range='/home/xfu/Project/CZQ/RNA-Seq/track/E2-1.bedgraph', name='E', col.mountain='white', fill.mountain=c('white', col_set[4]))
dtrack4.4 <- DataTrack(range='/home/xfu/Project/CZQ/RNA-Seq/track/D2-1.bedgraph', name='D', col.mountain='white', fill.mountain=c('white', col_set[4]))

save(list = c('txdb', 'grtrack', 'gtrack', 
              'dtrack1.1', 'dtrack1.2', 'dtrack1.3', 'dtrack1.4', 
              'dtrack2.1', 'dtrack2.2', 'dtrack2.3', 'dtrack2.4',
              'dtrack3.1', 'dtrack3.2', 'dtrack3.3', 'dtrack3.4',
              'dtrack4.1', 'dtrack4.2', 'dtrack4.3', 'dtrack4.4'),
              file='RData/tracks.RData')
