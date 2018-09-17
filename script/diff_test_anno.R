suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(DiffBind))
library(xlsx)

argv <- commandArgs(T)
input_rdata <- argv[1]
sqlite <- argv[2]
output_rdata <- argv[3]
output_xlsx <- argv[4]

load(input_rdata)
txdb <- loadDb(sqlite)

gene_anno <- function(x) {
  peakAnno <- annotatePeak(x, TxDb=txdb, tssRegion = c(-2000,0))
  anno <- as.data.frame(peakAnno)
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
