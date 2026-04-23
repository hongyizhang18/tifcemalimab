# This script visualizes immune infiltration from transcriptomic data.
# Updated: 2025-09-27

# Load packages
ps = c('tidyverse','openxlsx','clusterProfiler','ggplot2','circlize','pheatmap','ComplexHeatmap')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Immune cell data import
## Immune cell scores
immune = read.csv('immune.csv', header = TRUE, row.names = 1)
subtype = read.csv('ccsubtype.csv', header = TRUE)
estimate = immune %>% dplyr::filter(grepl('estimate', rownames(.))) %>%
  t() %>% scale() %>% as.data.frame() %>% rownames_to_column('sample_id')
## Clinical metadata
meta = read.xlsx('sample-info.xlsx', sheet = 1)
meta = meta %>% dplyr::filter(sample_type %in% c('T'))
table(colnames(immune) == meta$sample_id)

#### Bubble plot: pre vs post
# Reshape to long format
compare = immune %>% dplyr::filter(grepl('xcell', rownames(.))) %>%
  dplyr::filter(!grepl('score|value', rownames(.))) %>%
  rownames_to_column('obj') %>%
  pivot_longer(., -obj, names_to = 'sample_id', values_to = 'expr') %>%
  merge(., meta %>% dplyr::select(c('sample_id','sample_group')), by = 'sample_id')
colnames(compare) = c('sample_id','obj','expr','group')
compare$group = factor(compare$group, levels = c('combo.T.post','combo.T.pre','single.T.post','single.T.pre'))
# Differential analysis across groups
res = compare %>% group_by(obj) %>%
  summarise(combo_pval = wilcox.test(expr[group == levels(compare$group)[1]], expr[group == levels(compare$group)[2]], paired = FALSE, exact = FALSE)$p.value,
            single_pval = wilcox.test(expr[group == levels(compare$group)[3]], expr[group == levels(compare$group)[4]], paired = FALSE, exact = FALSE)$p.value,
            combo_logfc = log2(mean(expr[group == levels(compare$group)[1]]) / mean(expr[group == levels(compare$group)[2]])),
            single_logfc = log2(mean(expr[group == levels(compare$group)[3]]) / mean(expr[group == levels(compare$group)[4]]))) %>%
  mutate(combo_padj = round(p.adjust(combo_pval, method = 'BH'), 10),
         single_padj = round(p.adjust(single_pval, method = 'BH'), 10))
# Bubble plot for group differences
bubble = res %>% dplyr::select(c('obj','combo_logfc','combo_padj','single_logfc','single_padj')) %>%
  pivot_longer(cols = -obj, names_to = c('comp', '.value'), names_sep = '_') %>%
  mutate(obj = str_remove(obj, '^[^.]*\\.'),
         comp = factor(comp, levels = c('combo', 'single')),
         color = ifelse(logfc > 0, -log10(padj), -(-log10(padj))),
         size = abs(logfc)) %>%
  dplyr::filter(!is.infinite(logfc) & !is.na(logfc) & logfc != 'NaN')
summary((bubble$size))
summary((bubble$color))
ggplot(bubble, aes(x = obj, y = comp)) +
  geom_point(aes(size = size, color = color), alpha = 0.9) +
  geom_text(data = subset(bubble, padj < 0.05), aes(label = '*'), color = 'black', size = 5, vjust = 0.7) +
  scale_size_continuous(name = 'log2fc', range = c(2,8), breaks = seq(0, 5, length.out = 5)) +
  scale_color_gradient2(name = '-log10(adj.P)', low = '#0096c7', mid = 'white', high = '#d00000', midpoint = 0, limits = c(-4, 4)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = 'grey90', linewidth = 0.15),
        legend.position = 'none')
# Bubble plot for cross-group comparison
bubble = res %>% dplyr::select(c('obj','combo_logfc','combo_padj','single_logfc','single_padj')) %>%
  mutate(obj = str_remove(obj, '^[^.]*\\.'),
         sig = ifelse(combo_padj < 0.05 & single_padj < 0.05, 'both', ifelse(combo_padj < 0.05, 'combo only', ifelse(single_padj < 0.05, 'single only', 'not sig'))))
summary(bubble$combo_logfc)
summary(bubble$single_logfc)
ggplot(bubble, aes(x = combo_logfc, y = single_logfc)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dotted', color = 'gray60', alpha = 0.5) +
  geom_point(aes(color = sig), size = 3, alpha = 0.7) +
  geom_text_repel(data = subset(bubble, sig != 'not sig'), aes(label = obj), size = 3, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', segment.size = 0.5, segment.alpha = 0.6, force = 1) +
  scale_color_manual(name = 'significance', values = c('both' = '#E41A1C', 'combo only' = '#377EB8', 'single only' = '#4DAF4A', 'not sig' = 'gray70')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'grey90', linewidth = 0.15),
        legend.position = 'none') +
  xlim(-5, 5) + ylim(-5, 5) +
  labs(x = 'combo logFC', y = 'single logFC')

#### Bubble plot: responder vs non-responder
# Reshape to long format
compare = immune %>% dplyr::filter(grepl('xcell', rownames(.))) %>%
  dplyr::filter(!grepl('score|value|Correlation', rownames(.))) %>%
  rownames_to_column('obj') %>%
  pivot_longer(., -obj, names_to = 'sample_id', values_to = 'expr') %>%
  merge(., meta %>% dplyr::select(c('sample_id','sample_group')), by = 'sample_id')
colnames(compare) = c('sample_id','obj','expr','group')
compare$group = factor(compare$group, levels = c('combo.T.mpr','combo.T.nmpr','single.T.mpr','single.T.nmpr'))
compare$group = factor(compare$group, levels = c('combo.T.cpr','combo.T.ncpr','single.T.cpr','single.T.ncpr'))
# Differential analysis across groups
res = compare %>% filter(!is.na(group)) %>% group_by(obj) %>%
  summarise(combo_mw = wilcox.test(expr[group == levels(compare$group)[1]], expr[group == levels(compare$group)[2]], paired = FALSE, exact = FALSE)$p.value,
            single_mw = wilcox.test(expr[group == levels(compare$group)[3]], expr[group == levels(compare$group)[4]], paired = FALSE, exact = FALSE)$p.value,
            combo_ks = ks.test(expr[group == levels(compare$group)[1]], expr[group == levels(compare$group)[2]])$p.value,
            single_ks = ks.test(expr[group == levels(compare$group)[3]], expr[group == levels(compare$group)[4]])$p.value,
            combo_logfc = log2(mean(expr[group == levels(compare$group)[1]]) / mean(expr[group == levels(compare$group)[2]])),
            single_logfc = log2(mean(expr[group == levels(compare$group)[3]]) / mean(expr[group == levels(compare$group)[4]]))) %>%
  mutate(combo_mw_padj = round(p.adjust(combo_mw, method = 'BH'), 10),
         single_mw_padj = round(p.adjust(single_mw, method = 'BH'), 10),
         combo_ks_padj = round(p.adjust(combo_ks, method = 'BH'), 10),
         single_ks_padj = round(p.adjust(single_ks, method = 'BH'), 10))
# Bubble plot for group differences
bubble = res %>% dplyr::select(c('obj','combo_logfc','combo_mw','single_logfc','single_mw')) %>%
  pivot_longer(cols = -obj, names_to = c('comp', '.value'), names_sep = '_') %>%
  mutate(obj = str_remove(obj, '^[^.]*\\.'), comp = factor(comp, levels = c('combo', 'single')),
         color = ifelse(logfc > 0, -log10(mw), -(-log10(mw))), size = abs(logfc)) %>%
  dplyr::filter(!is.infinite(logfc) & !is.na(logfc) & logfc != 'NaN')
summary((bubble$size))
summary((bubble$color))
ggplot(bubble, aes(x = obj, y = comp)) +
  geom_point(aes(size = size, color = color), alpha = 0.9) +
  geom_text(data = subset(bubble, mw < 0.1), aes(label = '*'), color = 'black', size = 5, vjust = 0.7) +
  scale_size_continuous(name = 'log2fc', range = c(2,8), breaks = seq(0, 50, length.out = 5)) +
  scale_color_gradient2(name = '-log10(adj.P)', low = '#0096c7', mid = 'white', high = '#d00000', midpoint = 0, limits = c(-3, 3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = 'grey90', linewidth = 0.15),
        legend.position = 'none')
# Bubble plot for cross-group comparison
bubble = res %>% dplyr::select(c('obj','combo_logfc','combo_mw_padj','single_logfc','single_mw_padj')) %>%
  mutate(obj = str_remove(obj, '^[^.]*\\.'),
         sig = ifelse(combo_mw_padj < 0.05 & single_mw_padj < 0.05, 'both', ifelse(combo_mw_padj < 0.05, 'combo only', ifelse(single_mw_padj < 0.05, 'single only', 'not sig'))))
summary(bubble$combo_logfc)
summary(bubble$single_logfc)
ggplot(bubble, aes(x = combo_logfc, y = single_logfc)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dotted', color = 'gray60', alpha = 0.5) +
  geom_point(aes(color = sig), size = 3, alpha = 0.7) +
  geom_text_repel(data = subset(bubble, sig != 'not sig'), aes(label = obj), size = 3, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', segment.size = 0.5, segment.alpha = 0.6, force = 1) +
  scale_color_manual(name = 'significance', values = c('both' = '#E41A1C', 'combo only' = '#377EB8', 'single only' = '#4DAF4A', 'not sig' = 'gray70')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'grey90', linewidth = 0.15),
        legend.position = 'none') +
  xlim(-50, 50) + ylim(-50, 50) +
  labs(x = 'combo logFC', y = 'single logFC')


#### Immune cell scores浸润热图
# Integrate data
anno = meta %>% dplyr::select(c('sample_id','sample_group','recist','cpr')) %>% 
  merge(., subtype %>% dplyr::select(c('sample_id','immune')), by = 'sample_id')
colnames(anno) = c('sample_id','group','recist','cpr','immune')
anno$group = factor(anno$group, levels = c('combo.T.pre','combo.T.post','single.T.pre','single.T.post'))
anno$recist = factor(anno$recist, levels = c('CR','PR','SD','PD'))
# Visualize heatmap
anno = anno %>% arrange(group, recist, immune) %>% column_to_rownames('sample_id')
df = immune %>% dplyr::filter(grepl('xcell', rownames(.))) %>% dplyr::filter(!grepl('score|value', rownames(.)))
df = df %>% dplyr::select(rownames(anno))

# Cell infiltration heatmap
rowSums(df)
bk = c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
col = HeatmapAnnotation(
  df = anno[, c('group','recist','cpr','immune')],
  col = list(
    group = c('combo.T.pre'='#f2cc8f', 'combo.T.post'='#f07167', 'single.T.pre' = '#57cc99', 'single.T.post' = '#0081a7'),
    recist = c('CR'='#e29578', 'PR'='#ffdab9', 'SD'='#83c5be', 'PD'='#006d77'),
    cpr = c('CPR'='#00a5cf', 'non-CPR'='#accbe1'),
    immune = c('Immune-Enriched, Non-Fibrotic'='#bc4749', 'Immune-Enriched, Fibrotic'='#a7c957', 'Fibrotic'='#6a994e', 'Depleted'='#386641')),
  show_legend = TRUE,
  annotation_name_side = 'left')
ht = Heatmap(
  as.matrix(df), top_annotation = col,
  col = c(colorRampPalette(colors = c('#1e96fc','#f8f9fa'))(length(bk)/2),
          colorRampPalette(colors = c('#f8f9fa','#f77f00'))(length(bk)/2)),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = T, show_column_names = F, border = FALSE,
  row_names_gp = gpar(fontsize = 10), row_names_side = 'left')

# Bubble plot annotation
bubble = res %>% dplyr::select(c('obj','combo_logfc','single_logfc')) %>% column_to_rownames('obj')
bubble = bubble[rownames(df), , drop = FALSE]
row = Heatmap(
  as.matrix(bubble),
  col = colorRamp2(c(-2, 0, 2), c('#1982c4', '#FFFFFF', '#c44536')),
  rect_gp = gpar(type = 'none'),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(x = x, y = y, r = unit(0.1, 'npc'),
                gp = gpar(fill = fill, col = NA))},
  cluster_rows = FALSEALSE, cluster_columns = FALSEALSE,
  show_row_names = F, show_column_names = T, show_heatmap_legend = T,
  width = unit(1, 'cm'))
draw(ht + row, ht_gap = unit(0.2, 'cm'))

