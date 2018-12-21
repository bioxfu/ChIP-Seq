library(Gviz)
library(GenomicFeatures)
library(RColorBrewer)
col_set <- brewer.pal(9, 'Set1')
options(ucscChromosomeNames=FALSE)
options(stringsAsFactors=FALSE)

load('RData/tracks.RData')
genes <- read.table('/home/xfu/Gmatic6/genome/PhytozomeV12/Bdistachyon/annotation/Bdistachyon_314_v3.1.gene.bed', row.names=4)

# DNA methylation
id_lst <- c(
  'Bradi1g38325',
  'Bradi1g45760', 
  'Bradi3g33130', 
  'Bradi1g45852', 
  'Bradi3g22930', 
  'Bradi1g51554',
  'Bradi2g28421'
  )
for (id in id_lst) {
  pdf(paste0('figure/track/', id, '_MeDIP_tracks.pdf'), wid=4, hei=7)
  plotTracks(list(gtrack, dtrack1.1, dtrack1.2, dtrack1.3, dtrack1.4, 
                  dtrack4.1, dtrack4.2, dtrack4.3, dtrack4.4, grtrack),
             chromosome=genes[id, 1], from=genes[id, 2] - 2000, to=genes[id, 3] + 2000, lwd=2, type='polygon',
             background.title='white', fontcolor.title='black', col.axis='black', cex.title=1, cex.axis=0.7, cex=1, cex.group=1)
  dev.off()
}

# H3K4me3
id_lst <- c(
  'Bradi1g53654',
  'Bradi2g33682',
  'Bradi2g59119'
)
for (id in id_lst) {
  pdf(paste0('figure/track/', id, '_H3K4me3_tracks.pdf'), wid=4, hei=7)
  plotTracks(list(gtrack, dtrack2.1, dtrack2.2, dtrack2.3, dtrack2.4, 
                  dtrack4.1, dtrack4.2, dtrack4.3, dtrack4.4, grtrack),
             chromosome=genes[id, 1], from=genes[id, 2] - 2000, to=genes[id, 3] + 2000, lwd=2, type='polygon',
             background.title='white', fontcolor.title='black', col.axis='black', cex.title=1, cex.axis=0.7, cex=1, cex.group=1)
  dev.off()
}

# H3K9me3
id_lst <- c(
  'Bradi3g59147',
  'Bradi3g31777'
)
for (id in id_lst) {
  pdf(paste0('figure/track/', id, '_H3K9me3_tracks.pdf'), wid=4, hei=7)
  plotTracks(list(gtrack, dtrack3.1, dtrack3.2, dtrack3.3, dtrack3.4, 
                  dtrack4.1, dtrack4.2, dtrack4.3, dtrack4.4, grtrack),
             chromosome=genes[id, 1], from=genes[id, 2] - 2000, to=genes[id, 3] + 2000, lwd=2, type='polygon',
             background.title='white', fontcolor.title='black', col.axis='black', cex.title=1, cex.axis=0.7, cex=1, cex.group=1)
  dev.off()
}
