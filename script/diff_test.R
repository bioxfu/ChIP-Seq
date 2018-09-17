suppressPackageStartupMessages(library(DiffBind))

argv <- commandArgs(T)
input_rdata <- argv[1]
MASK <- argv[2]
output_rdata <- argv[3]
load(input_rdata)

# MASK <- 'H3K9me3,H3K4me3:A,B,C,D,E'
# MASK <- 'WT,lec1l1l'

mask_lst <- unlist(strsplit(MASK, '^', fixed = TRUE))
diff_peak <- list()
diff_name <- c()

for (m in mask_lst) {
  maskStr <- unlist(strsplit(m, ':'))
  factors <- strsplit(maskStr[1], ',')[[1]]

  for (j in 1:(length(factors)-1)) {
    for (k in (j+1):length(factors)) {
      g1 <- factors[j]
      g2 <- factors[k]
      cat(g1, g2, '\n')
      # mask
      g1_mask <- d_count$masks[[g1]]
      g2_mask <- d_count$masks[[g2]]
      print(g1_mask)
      print(g2_mask)
      # Establishing a contrast
      d_count2 <- dba.contrast(d_count, g1_mask, g2_mask, g1, g2)
      # Performing the differential analysis
      d_count_test <- dba.analyze(d_count2, method=DBA_DESEQ2, bCorPlot = F)
      # Retrieving the differentially bound sites
      d_report <- dba.report(d_count_test, method=DBA_DESEQ2, th=0.05, bCalled=T, bUsePval=TRUE)
      diff_peak <- c(diff_peak, d_report)
      diff_name <- c(diff_name, paste0(g1, '_', g2))
    }    
  }
  names(diff_peak) <- diff_name
}

save(diff_peak, file = output_rdata)

