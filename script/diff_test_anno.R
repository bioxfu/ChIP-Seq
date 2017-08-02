suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(DiffBind))
library(xlsx)

argv <- commandArgs(T)
input_rdata <- argv[1]
sqlite <- argv[2]
gene2GO <- argv[3]
output_rdata <- argv[4]
output_xlsx <- argv[5]

load(input_rdata)
txdb <- loadDb(sqlite)
GO <- read.table(gene2GO, sep='\t', quote="", stringsAsFactors = F)
colnames(GO) <- c('geneId', 'GOterm')

gene_anno <- function(x) {
  peakAnno <- annotatePeak(x, TxDb=txdb)
  anno <- as.data.frame(peakAnno)
  anno <- merge(anno, GO, all.x=T)
  return(anno)
}

diff_peak_anno <- lapply(diff_peak, gene_anno)
save(diff_peak_anno, file = output_rdata)

wb <- createWorkbook()
for (i in 1:length(diff_peak_anno)) {
  sheet <- createSheet(wb, sheetName=names(diff_peak_anno)[i])
  addDataFrame(diff_peak_anno[[i]], sheet, row.names=FALSE)
}
saveWorkbook(wb, output_xlsx)
