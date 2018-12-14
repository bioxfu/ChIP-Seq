suppressPackageStartupMessages(library(DiffBind))

argv <- commandArgs(T)
sheet <- argv[1]
# heatmap_peak <- argv[2]
# heatmap_count <- argv[3]
# PCA <- argv[4]
RData <- argv[2]

# reading in the peaksets
d_peak <- dba(sampleSheet = sheet, bCorPlot = F)

# correlation heatmap using occupance (peak caller score) data
# pdf(heatmap_peak)
# plot(d_peak)
# dev.off()

# counting reads
d_count <- dba.count(d_peak, bCorPlot = F)

# correlation heatmap using affinity (read count) data
# pdf(heatmap_count)
# plot(d_count)
# dev.off()

## PCA plot
# pdf(PCA)
# dba.plotPCA(d_count, DBA_FACTOR, label = DBA_REPLICATE)
# dev.off()

save(list = c('d_peak', 'd_count'), file = RData)

