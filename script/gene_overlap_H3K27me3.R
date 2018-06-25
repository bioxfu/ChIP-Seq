library(magrittr)
library(GeneOverlap)
library(gplots)
library(plotrix)

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

options(stringsAsFactors = F)
genome <- read.table('/home/xfu/Gmatic5/old_version/Gmatic3/bin/geneOverlap_bg_tair10_genes')
targets <- read.table(input, head=T, sep='\t', quote='')
targets <- targets[(targets$distanceToTSS < targets$geneLength & targets$distanceToTSS > -1000), ]

#k27_cxf <- read.table('/home/xfu/Gmatic5/pub_data/arab_H3K27me3/H3K27me3_CXF.txt')
#colnames(k27_cxf) <- 'geneId'

k27_zyj <- read.table('/home/xfu/Gmatic5/pub_data/arab_H3K27me3/H3K27me3_ZYJ.txt')
colnames(k27_zyj) <- 'geneId'

gene_lst <- list(targets=targets$geneId %>% unique(),
                 H3K27me3.ZYJ=k27_zyj$geneId %>% unique())

mat <- NULL
for (i in 1:(length(gene_lst)-1)) {
  for (j in (i+1):length(gene_lst)) {
    go.obj <- newGeneOverlap(gene_lst[[i]], gene_lst[[j]], genome.size=nrow(genome)) %>% testGeneOverlap()
    num <- getIntersection(go.obj) %>% length()
    p <- getPval(go.obj) %>% format(digits=2)
    mat <- rbind(mat, c(names(gene_lst)[i], names(gene_lst)[j], num, p))
  }
}
colnames(mat) <- c('List1', 'List2', 'Number', 'P-value')

pdf(output, wid=10)
layout(matrix(c(1,2),nrow=1), wid=c(2,1))
par(mar=c(2,2,2,2))
par(xpd=T)
venn(gene_lst)
plot(0, type='n', xaxt='n', yaxt='n', bty='n')
addtable2plot(0.45, 0, mat, display.rownames=F, display.colnames=T, bty='n')
dev.off()

