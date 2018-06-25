library(RColorBrewer)
col_set <- brewer.pal(4, 'Set1')

argv <- commandArgs(T)
chip_sample <- argv[1]
input_sample <- argv[2]
output <- argv[3]

chip <- as.matrix(read.table(chip_sample, sep='\t'))
chip <- chip[, c(-1,-2)]
lab <- c(which(chip[1,]=='-1.0Kb'), which(chip[1,]=='TSS'), which(chip[1,]=='TES'), which(chip[1,]=='1.0Kb'))
input <- as.matrix(read.table(input_sample, sep='\t'))
input <- input[, c(-1,-2)]

pdf(output, wid=5, hei=4)
ymin <- 2
ymax <- 3.5
plot(chip[2,], chip[3,], type='n', ylim=c(ymin, ymax), xlab='', ylab='Normalized read coverage', xaxt='n', yaxt='n', cex.lab=1.2)
abline(h=seq(ymin, ymax, 0.5), col='gray', lwd=0.5)
abline(v=seq(0, 300, 50), col='gray', lwd=0.5)
lines(chip[2,], chip[3,], lwd=3, col=col_set[1])
lines(input[2,], input[3,], lwd=3, col=col_set[2])
axis(1, lab, c('-1000', 'start', 'end', '1000'), cex.axis=1.2)
axis(2, seq(ymin, ymax, 0.5), seq(ymin, ymax, 0.5), las=2, cex.axis=1.2)
legend('topleft', c('ChIP', 'input'), col=col_set[1:2], bty='n', lwd=3, cex=1.2)
dev.off()
