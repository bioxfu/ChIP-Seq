# setting the java.parameters option before the rJava package is loaded.
options(java.parameters = "-Xmx8000m")
library(yaml)
library(xlsx)

argv <- commandArgs(T)
config_file <- argv[1]
input_dir <- argv[2]
output_xlsx <- argv[3]

config <- yaml.load_file(config_file)
files <- paste0(input_dir, '/', config$peak1, '_vs_', config$peak2, '_', config$p , '/_all_peak_MAvalues.anno')

wb <- createWorkbook()
for (i in 1:length(files)) {
  sheet <- createSheet(wb, sheetName=paste0(config$peak1, '_vs_', config$peak2)[i])
  dfm <- read.table(files[i], sep='\t', header=T, quote = "")
  dfm <- dfm[dfm$P_value < 0.05,]
  addDataFrame(dfm, sheet, row.names=FALSE)
}
saveWorkbook(wb, output_xlsx)
