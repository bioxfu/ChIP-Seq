peak <- 'peak_1e-2'

depth <- read.table('stat/bamqc_stat.tsv', header = T, sep = '\t')[4:5, 3]
coverage <- read.table(paste0(peak, '/peaks_olp_count_coverage'))
#coverage <- coverage[coverage$V4==2,]
IP_1 <- log2(coverage$V6 / depth[1] * 1000000)
IP_2 <- log2(coverage$V7 / depth[2] * 1000000)

r1 <- round(cor(IP_1, IP_2), 2)

pdf(paste0('figure/', peak, '_rep_correlations_PCC.pdf'))
plot(IP_1, IP_2, pch=16, cex=0.5, xlab='IP replicate 1', ylab='IP replicate 2', cex.axis=1.5, cex.lab=1.5)
text(5, 11, bquote(italic('r')==.(r1)), cex=1.5)
dev.off()
