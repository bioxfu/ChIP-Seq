suppressPackageStartupMessages(library(DiffBind))

argv <- commandArgs(T)
input_rdata <- argv[1]
output_pdf <- argv[2]
bedout <- argv[3]
submitout <- argv[4]

load(input_rdata)

mask <- unlist(strsplit(output_pdf, '_'))[3]

grange2bed <- function(x) {
  df <- data.frame(seqnames=seqnames(x),
                   starts=start(x)-1,
                   ends=end(x),
                   names=paste0('MACS_peak_', 1:length(x)),
                   scores=rep(".", length(x)))
  return(df)
}

grange2submit <- function(x) {
  df <- data.frame(seqnames=seqnames(x),
                   starts=start(x) + round((end(x) - start(x) + 1)/2) - 1,
                   ends=start(x) + round((end(x) - start(x) + 1)/2),
                   names=paste0('MACS_submit_', 1:length(x)),
                   scores=rep(".", length(x)))
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
  } else {
    print('Too many peaksets (n>4)')
  }
}

pdf(output_pdf)
find_overlap(mask)
dev.off()

