library(Gviz)
library(GenomicFeatures)
library(RColorBrewer)
col_set <- brewer.pal(9, 'Set1')
options(ucscChromosomeNames=FALSE)
options(stringsAsFactors=FALSE)

load('RData/tracks.RData')
genes <- read.table('/home/xfu/Gmatic6/genome/PhytozomeV12/Bdistachyon/annotation/Bdistachyon_314_v3.1.gene.bed', row.names=4)


# DNA methylation
track1_height_lim <- c(5, 2, 4, 2, 2, 4, 2)
track4_height_lim <- c(2, 25, 2, 0.5, 5, 3, 1.5)

id_lst <- c(
  'Bradi1g38325',
  'Bradi1g45760', 
  'Bradi3g33130', 
  'Bradi1g45852', 
  'Bradi3g22930', 
  'Bradi1g51554',
  'Bradi2g28421'
)
for (i in 1:length(id_lst)) {
  id <- id_lst[i]
  displayPars(dtrack1.1) <- list(ylim=c(0,track1_height_lim[i]))
  displayPars(dtrack1.2) <- list(ylim=c(0,track1_height_lim[i]))
  displayPars(dtrack1.3) <- list(ylim=c(0,track1_height_lim[i]))
  displayPars(dtrack1.4) <- list(ylim=c(0,track1_height_lim[i]))
  displayPars(dtrack4.1) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.2) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.3) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.4) <- list(ylim=c(0,track4_height_lim[i]))
  
  pdf(paste0('figure/track/final/', id, '_MeDIP_tracks_final.pdf'), wid=4, hei=7)
  plotTracks(list(gtrack, dtrack1.1, dtrack1.2, dtrack1.3, dtrack1.4, 
                  dtrack4.1, dtrack4.2, dtrack4.3, dtrack4.4, grtrack),
             chromosome=genes[id, 1], from=genes[id, 2] - 2000, to=genes[id, 3] + 2000, lwd=2, type='polygon',
             background.title='white', fontcolor.title='black', col.axis='black', cex.title=1, cex.axis=0.7, cex=1, cex.group=1)
  dev.off()
}

# H3K4me3
track2_height_lim <- c(2, 1.5, 1)
track4_height_lim <- c(3, 50, 2)

id_lst <- c(
  'Bradi1g53654',
  'Bradi2g33682',
  'Bradi2g59119'
)
for (i in 1:length(id_lst)) {
  id <- id_lst[i]
  displayPars(dtrack2.1) <- list(ylim=c(0,track2_height_lim[i]))
  displayPars(dtrack2.2) <- list(ylim=c(0,track2_height_lim[i]))
  displayPars(dtrack2.3) <- list(ylim=c(0,track2_height_lim[i]))
  displayPars(dtrack2.4) <- list(ylim=c(0,track2_height_lim[i]))
  displayPars(dtrack4.1) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.2) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.3) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.4) <- list(ylim=c(0,track4_height_lim[i]))

  pdf(paste0('figure/track/final/', id, '_H3K4me3_tracks_final.pdf'), wid=4, hei=7)
  plotTracks(list(gtrack, dtrack2.1, dtrack2.2, dtrack2.3, dtrack2.4, 
                  dtrack4.1, dtrack4.2, dtrack4.3, dtrack4.4, grtrack),
             chromosome=genes[id, 1], from=genes[id, 2] - 2000, to=genes[id, 3] + 2000, lwd=2, type='polygon',
             background.title='white', fontcolor.title='black', col.axis='black', cex.title=1, cex.axis=0.7, cex=1, cex.group=1)
  dev.off()
}


# H3K9me3
track3_height_lim <- c(1, 1.5)
track4_height_lim <- c(1, 50)

id_lst <- c(
  'Bradi3g59147',
  'Bradi3g31777'
)

for (i in 1:length(id_lst)) {
  id <- id_lst[i]
  displayPars(dtrack3.1) <- list(ylim=c(0,track3_height_lim[i]))
  displayPars(dtrack3.2) <- list(ylim=c(0,track3_height_lim[i]))
  displayPars(dtrack3.3) <- list(ylim=c(0,track3_height_lim[i]))
  displayPars(dtrack3.4) <- list(ylim=c(0,track3_height_lim[i]))
  displayPars(dtrack4.1) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.2) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.3) <- list(ylim=c(0,track4_height_lim[i]))
  displayPars(dtrack4.4) <- list(ylim=c(0,track4_height_lim[i]))

  pdf(paste0('figure/track/final/', id, '_H3K9me3_tracks_final.pdf'), wid=4, hei=7)
  plotTracks(list(gtrack, dtrack3.1, dtrack3.2, dtrack3.3, dtrack3.4, 
                  dtrack4.1, dtrack4.2, dtrack4.3, dtrack4.4, grtrack),
             chromosome=genes[id, 1], from=genes[id, 2] - 2000, to=genes[id, 3] + 2000, lwd=2, type='polygon',
             background.title='white', fontcolor.title='black', col.axis='black', cex.title=1, cex.axis=0.7, cex=1, cex.group=1)
  dev.off()
}
