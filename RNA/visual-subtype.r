# This script visualizes transcriptomic subtype clustering.
# Updated: 2026-01-17

# Load packages
ps = c('tidyverse','GSVA','GSEABase','clusterProfiler','pheatmap','ggplot2','ggalluvial')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Data preprocessing
# Clinical metadata
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>% dplyr::filter(sample_type %in% c('T'))
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1', immune == 'Immune-Enriched, Fibrotic' ~ 'C2', immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'))
meta = merge(meta, subtype, by = 'sample_id')
meta$treatment = factor(meta$treatment, levels = c('combo','single'))
meta$immune = factor(meta$immune, levels = c('C1','C2','C3','C4'))
meta$recist = factor(meta$recist, levels = c('CR','PR','SD','PD'))
meta$response = factor(meta$response, levels = c('cpr','mpr','non-mpr','non-surgery'))
anno = meta %>% arrange(immune, treatment, recist, response) %>% dplyr::select(sample_id, immune, treatment, recist, response) %>% column_to_rownames('sample_id')

# Transcriptome matrix
matrix = read.csv('tpm.csv', header = TRUE, row.names = 1)
matrix = log2(matrix + 1)
# Gene sets
sig.tme = getGmt(file.path(data_dir, '29TME.gmt'), geneIdType = SymbolIdentifier())
sig.hmk = getGmt(file.path(data_dir, 'h.all.2024.Hs.gmt'), geneIdType = SymbolIdentifier())
# Signature scoring
scores.tme = gsva(as.matrix(matrix), sig.tme, method = 'ssgsea', kcdf = 'Gaussian', min.sz = 1, max.sz = Inf) %>% as.data.frame() %>% dplyr::select(rownames(anno))
scores.hmk = gsva(as.matrix(matrix), sig.hmk, method = 'ssgsea', kcdf = 'Gaussian', min.sz = 1, max.sz = Inf) %>% as.data.frame() %>% dplyr::select(rownames(anno))
rownames(scores.hmk) = gsub('HALLMARK_', '', rownames(scores.hmk))
scores.hmk = scores.hmk[c('ADIPOGENESIS','CHOLESTEROL_HOMEOSTASIS','FATTY_ACID_METABOLISM','GLYCOLYSIS','OXIDATIVE_PHOSPHORYLATION',
                          'E2F_TARGETS','G2M_CHECKPOINT','MYC_TARGETS_V2','DNA_REPAIR','ANGIOGENESIS','EPITHELIAL_MESENCHYMAL_TRANSITION','HYPOXIA',
                          'ALLOGRAFT_REJECTION','INFLAMMATORY_RESPONSE','INTERFERON_ALPHA_RESPONSE','INTERFERON_GAMMA_RESPONSE','IL2_STAT5_SIGNALING','TNFA_SIGNALING_VIA_NFKB','TGF_BETA_SIGNALING',
                          'NOTCH_SIGNALING','HEDGEHOG_SIGNALING','PI3K_AKT_MTOR_SIGNALING','MTORC1_SIGNALING','KRAS_SIGNALING_UP','KRAS_SIGNALING_DN'), ]

#### Data visualization
# Main figure
df = scores.hmk %>% t() %>% as.data.frame() %>% scale() %>% t() %>% as.data.frame()
df = scores.hmk %>% apply(., 1, function(x) (x - mean(x)) / sd(x)) %>% t()
rowSums(df)
bk = c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(df, annotation_col = anno, 
         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
         color = c(colorRampPalette(colors = c('#023e8a','#8ecae6'))(length(bk)/4),
                   colorRampPalette(colors = c('#8ecae6','#f8f9fa'))(length(bk)/4),
                   colorRampPalette(colors = c('#f8f9fa','#f4a261'))(length(bk)/4),
                   colorRampPalette(colors = c('#f4a261','#ba5624'))(length(bk)/4)),
         annotation_colors = list(
           immune = c('C1'='#bc4749', 'C2'='#eec170', 'C3'='#a7c957', 'C4'='#0D6E6E'),
           treatment = c('combo'='#a978a3', 'single'='#f4c0cd'),
           recist = c('CR'='#caf0f8', 'PR'='#90e0ef', 'SD'='#00b4d8', 'PD'='#0077b6'),
           response = c('cpr'='#ed6a5a', 'mpr'='#ffd166', 'non-mpr'='#5ca4a9', 'non-surgery'='#adb5bd')))
# Supplementary figure
df = scores.tme %>% t() %>% as.data.frame() %>% scale() %>% t() %>% as.data.frame()
df = scores.tme %>% apply(., 1, function(x) (x - mean(x)) / sd(x)) %>% t()
rowSums(df)
bk = c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(df, annotation_col = anno, 
         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
         color = c(colorRampPalette(colors = c('#023e8a','#8ecae6'))(length(bk)/4),
                   colorRampPalette(colors = c('#8ecae6','#f8f9fa'))(length(bk)/4),
                   colorRampPalette(colors = c('#f8f9fa','#f4a261'))(length(bk)/4),
                   colorRampPalette(colors = c('#f4a261','#ba5624'))(length(bk)/4)),
         annotation_colors = list(
           immune = c('C1'='#bc4749', 'C2'='#eec170', 'C3'='#a7c957', 'C4'='#0D6E6E'),
           treatment = c('combo'='#a978a3', 'single'='#f4c0cd'),
           recist = c('CR'='#caf0f8', 'PR'='#90e0ef', 'SD'='#00b4d8', 'PD'='#0077b6'),
           response = c('cpr'='#ed6a5a', 'mpr'='#ffd166', 'non-mpr'='#5ca4a9', 'non-surgery'='#adb5bd')))

