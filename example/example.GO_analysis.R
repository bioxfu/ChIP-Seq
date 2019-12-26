library(pheatmap)
library(RColorBrewer)
source('script/topGO.R')

folder <- 'peak_1e-2'

input <- dir(folder, pattern='_anno.xls')

go_lst <- list()

for (i in 1:length(input)) {

  dfm <- read.table(paste0(folder, '/', input[i]), header = T, quote = '', comment.char = '', sep = '\t', stringsAsFactors = F)

  gene <- unique(dfm$geneId[dfm$distanceToTSS >= -2000 & dfm$distanceToTSS <= dfm$geneLength])

  go <- topGO(gene, species = 'tair10')

  go_lst[[i]] <- go

}

names(go_lst) <- input

write2xlsx(go_lst, paste0(folder, '/target_gene_GO_enrichment.xlsx'))


