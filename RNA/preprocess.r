# This script converts gene IDs and removes duplicate genes.
# Updated: 2025-12-04

# Load packages
ps = c('tidyverse','data.table','clusterProfiler','org.Hs.eg.db','biomaRt','ggplot2','ggpubr')
for (i in ps) {
  library(i, character.only = TRUE)
}

#### Data import
## Expression matrix
df = read.table('rna_features.txt', header = TRUERUE)
meta = df[ ,c(1,6)]
expr = df[ ,7:ncol(df)]
colnames(meta) = c('ENSEMBL','Length')
colnames(expr) = sapply(strsplit(colnames(expr), '\\.'), `[`, 3)
expr = cbind(ENSEMBL = sub('\\..*$', '', meta$ENSEMBL), Length = meta$Length, expr)
table(duplicated(expr$ENSEMBL))

#### Gene ID conversion and deduplication
## Query biomaRt
mart = useMart('ensembl', dataset = 'hsapiens_gene_ensembl', host = 'https://asia.ensembl.org')
symbol = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = expr$ENSEMBL, mart = mart)
expr = merge(expr, symbol, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
expr = expr %>% filter(!is.na(external_gene_name) & external_gene_name != '')
table(duplicated(expr$external_gene_name))
table(is.na(expr$external_gene_name))
## Keep maximum expression per gene
expr.max = aggregate(. ~ external_gene_name, data = expr, FUN = max)
meta = expr.max[ ,1:3]
count = expr.max[ ,-c(2:3)]
count = column_to_rownames(count, var = 'external_gene_name')
## Save outputs
write.csv(count, file = 'count.csv', row.names = TRUE)
write.csv(meta, file = 'count-meta.csv', row.names = FALSE)

