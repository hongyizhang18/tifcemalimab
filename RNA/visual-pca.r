# This script visualizes principal component analysis for the expression matrix.
# Updated: 2025-12-04

# Load packages
ps = c('tidyverse','openxlsx','data.table','edgeR','FactoMineR','factoextra','ggplot2','ggrepel','ggpubr')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Data import
# Clinical metadata
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>% dplyr::filter(sample_type == 'T')
# Transcriptome matrix
matrix = read.csv('tpm.csv', header = TRUE, row.names = 1)

#### Data preprocessing
meta = meta %>% dplyr::filter(treatment %in% c('single'))
matrix = matrix %>% dplyr::select(meta$sample_id)

#### Principal component analysis
metric = meta$group
# Visualize with factoextra
data = matrix %>% t() %>% as.data.frame() %>% na.omit()
data$group = factor(metric, levels = unique(metric))
pca.data = PCA(data[, -ncol(data)], graph = FALSE)
palette = if (length(unique(metric)) > 4) c('#0d47a1','#1F78B4','#70b77e','#c5d86d','#e3bd64','#f4a261','#db6d6b','#e07ba7','#5e548e','#e0b1cb','#717596','#424b54') else c('#1F78B4','#70b77e','#e3bd64','#db6d6b')
fviz_pca_ind(pca.data, col.ind = data$group, 
             geom.ind = 'point',
             palette = palette, addEllipses = T) + theme_bw()

