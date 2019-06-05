suppressPackageStartupMessages(library(DiffBind))

argv <- commandArgs(T)
input_rdata <- argv[1]
output_pdf <- argv[2]
bedout <- argv[3]
submitout <- argv[4]
bedout2 <- argv[5]
submitout2 <- argv[6]

# input_rdata <- 'RData/peak_1e-5_dba.RData'
# output_pdf <- 'figure/peak_1e-5_lec1l1l:2.3_overlap.pdf'
# bedout <- 'peak_1e-5/lec1l1l:2.3_consensusPeaks.bed'
# submitout <- 'peak_1e-5/lec1l1l:2.3_consensusSubmits.bed'
# bedout2 <- 'peak_1e-5/lec1l1l:2.3_any2Peaks.bed'
# submitout2 <- 'peak_1e-5/lec1l1l:2.3_any2Submits.bed'
 
load(input_rdata)

mask <- unlist(strsplit(output_pdf, '_'))[3]

grange2bed <- function(x) {
  df <- data.frame(seqnames=seqnames(x),
                   #starts=start(x)-1,
                   starts=start(x),
                   ends=end(x),
                   names=paste0('MACS_peak_', 1:length(x)),
                   scores=apply(as.data.frame(x)[-c(1:5)], 1, max))
  # the scores associated with each site are derived from teh peak caller confidence score
  # and are measure of confidence in the peak call, not a measure of how strong or distince
  # the peak is
  return(df)
}

grange2submit <- function(x) {
  df <- data.frame(seqnames=seqnames(x),
                   #starts=start(x) + round((end(x) - start(x) + 1)/2) - 1,
                   starts=start(x) + round((end(x) - start(x) + 1)/2),
                   ends=start(x) + round((end(x) - start(x) + 1)/2),
                   names=paste0('MACS_submit_', 1:length(x)),
                   scores=apply(as.data.frame(x)[-c(1:5)], 1, max))
  return(df)
}

find_overlap <- function(x) {
  print(x)
  maskStr <- unlist(strsplit(x, ':'))
  factors <- strsplit(maskStr[1], ',')[[1]]
  print(factors)

  mask.factors <- d_peak$masks[[factors[1]]]
  for (i in factors) {
    mask.factors <- mask.factors | d_peak$masks[[i]]
  }

  mask <- mask.factors
  
  if (sum(mask) <= 4) {
    olp <- dba.overlap(d_peak, mask, mode=DBA_OLAP_PEAKS)
    olp_num <- unlist(lapply(olp, length))
    inAll_prop <- round(olp_num['inAll'] / sum(olp_num) * 100, 1)
    dba.plotVenn(d_peak, mask, main=x, sub=paste0('inAll: ', inAll_prop, '%'))
    write.table(grange2bed(olp$inAll), bedout, col.names = F, row.names = F, sep = '\t', quote = F)
    write.table(grange2submit(olp$inAll), submitout, col.names = F, row.names = F, sep = '\t', quote = F)
    if (sum(mask) == 3) {
      any2bed <- rbind(grange2bed(olp$inAll), grange2bed(olp$notA), grange2bed(olp$notB), grange2bed(olp$notC))
      any2submit <- rbind(grange2submit(olp$inAll), grange2submit(olp$notA), grange2submit(olp$notB), grange2submit(olp$notC))
    }
    else if (sum(mask) == 2) {
      any2bed <- grange2bed(olp$inAll)
      any2submit <- grange2submit(olp$inAll)
    }
    write.table(any2bed, bedout2, col.names = F, row.names = F, sep = '\t', quote = F)
    write.table(any2submit, submitout2, col.names = F, row.names = F, sep = '\t', quote = F)
  } else {
    print('Too many peaksets (n>4)')
  }
}

pdf(output_pdf)
find_overlap(mask)
dev.off()

