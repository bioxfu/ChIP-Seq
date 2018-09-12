library(rvest)

parse_html <- function(html) {
  html <- read_html(html)
  x <- html %>% html_nodes('.table-summary') %>% html_text()
  strsplit(x[3], '\n')[[1]][c(6,8)]
}

dfm <- data.frame(sample_name = sub('.bamqc', '', dir('bam', pattern = '*.bamqc')),
                  file = paste0(dir('bam', pattern = '*.bamqc', full=TRUE), '/qualimapReport.html'))

clean_reads <- unique(read.table('stat/fastqc_stat.tsv', head=T, row.names = 1)$tot.clean.seq) * 2

stat <- apply(dfm, 1, function(x){parse_html(x[2])}) %>% t()
stat <- as.data.frame(cbind(clean_reads, as.numeric(sub(' / 100%', '', gsub(',', '', stat[, 2])))))
rownames(stat) <- dfm$sample_name
colnames(stat) <- c('Number of clean reads', 'Unique mapped reads')
stat$Percent <- round(stat[, 2] / stat[, 1] * 100, 2)
write.table(stat, 'stat/bamqc_stat.tsv', quote=F, col.names = NA, sep='\t')
