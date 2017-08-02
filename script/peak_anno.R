suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ChIPseeker))

argv <- commandArgs(T)
input_bed <- argv[1]
sqlite <- argv[2]
gene2GO <- argv[3]
output_xls <- argv[4]
output_pie <- argv[5]

peak <- readPeakFile(input_bed)
txdb <- loadDb(sqlite)
GO <- read.table(gene2GO, sep='\t', quote="", stringsAsFactors = F)
colnames(GO) <- c('geneId', 'GOterm')

gene_anno <- function(x) {
  peakAnno <- annotatePeak(x, TxDb=txdb)
  anno <- as.data.frame(peakAnno)
  anno <- merge(anno, GO, all.x=T)
  return(list(peakAnno=peakAnno, anno=anno))
}

peak_anno <- gene_anno(peak) 
write.table(peak_anno$anno, output_xls, sep='\t', quote=F, row.names=F)

pdf(output_pie, width=7, height=5)
plotAnnoPie(peak_anno$peakAnno)
dev.off()
