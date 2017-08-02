suppressPackageStartupMessages(library(DiffBind))

argv <- commandArgs(T)
input_rdata <- argv[1]
MASK <- argv[2]
output_pdf <- argv[3]
load(input_rdata)

# MASK <- 'H3K4me3:A,B,C,D^H3K9me3:A,B,C,E'

mask_lst <- unlist(strsplit(MASK, '^', fixed = TRUE))

find_overlap <- function(x) {
  print(x)
  maskStr <- unlist(strsplit(x, ':'))
  factors <- strsplit(maskStr[1], ',')[[1]]
  treats <- strsplit(maskStr[2], ',')[[1]]
  
  mask.factors <- d_peak$masks[[factors[1]]]
  for (i in factors) {
    mask.factors <- mask.factors | d_peak$masks[[i]]
  }
  
  mask.treats <- d_peak$masks[[treats[1]]]
  for (i in treats) {
    mask.treats <- mask.treats | d_peak$masks[[i]]
  }
  
  mask <- mask.factors & mask.treats
  
  if (sum(mask) <= 4) {
    olp <- dba.overlap(d_peak, mask, mode=DBA_OLAP_PEAKS)
    olp_num <- unlist(lapply(olp, length))
    inAll_prop <- round(olp_num['inAll'] / sum(olp_num) * 100, 1)
    dba.plotVenn(d_peak, mask, main=x, sub=paste0('inAll: ', inAll_prop, '%'))
  } else {
    print('Too many peaksets (n>4)')
  }
}

pdf(output_pdf)
for (m in mask_lst) {
  find_overlap(m)
}
dev.off()

