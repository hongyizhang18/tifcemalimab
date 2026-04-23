# This script converts counts to TPM and evaluates batch-related variation.
# Updated: 2025-12-04

# Load packages
ps = c('tidyverse','openxlsx','sva','immunedeconv')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Data import
# Expression matrix
matrix.raw = read.csv('count.csv', header = TRUE, row.names = 1)
meta = read.xlsx('sample-info.xlsx', sheet = 1)
# Align sample identifiers
matrix.raw = matrix.raw %>% dplyr::select(meta$sample_name)
table(colnames(matrix.raw) == meta$sample_name)
colnames(matrix.raw) = meta$sample_id[match(colnames(matrix.raw), meta$sample_name)]
table(colnames(matrix.raw) == meta$sample_id)

#### Preprocessing
## TPM conversion
gene.length = read.csv('count-meta.csv', header = TRUE, row.names = 1)
table(rownames(matrix.raw) == rownames(gene.length))
str(matrix.raw)
count2tpm = function(counts, effLen) {
  rate = log(counts) - log(effLen)
  norm = log(sum(exp(rate)))
  tpm = exp(rate - norm + log(1e6))
  return(tpm)
}
table(is.na(gene.length$Length) | gene.length$Length == '')
matrix.tpm = as.data.frame(apply(matrix.raw, 2, count2tpm, effLen = as.numeric(gene.length$Length)))

#### Save expression matrices
# Count matrix
write.csv(matrix.raw, file = file.path(output_dir, 'count.csv'))
# TPM matrix
write.csv(matrix.tpm, file = file.path(output_dir, 'tpm.csv'))

# Quality control
estimate = deconvolute_estimate(matrix.tpm)
estimate = estimate %>% t() %>% scale() %>% as.data.frame() %>% rownames_to_column('sample_id') %>% dplyr::select(c('sample_id','TumorPurity'))
estimate = merge(estimate, meta, by = 'sample_id')
ggplot(estimate, aes(x = sample_type, y = TumorPurity, fill = sample_type)) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  theme_bw()

