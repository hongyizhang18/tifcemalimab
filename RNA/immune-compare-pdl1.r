# This script performs gene set scoring for transcriptomic data.
# Updated: 2025-09-27

# Load packages
ps = c('tidyverse','openxlsx','clusterProfiler','GSVA','ConsensusClusterPlus','pheatmap','ggplot2','ggrepel','ggpubr','ggalluvial')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Subtype boxplots
# Prepare data
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>%
  dplyr::filter(sample_type %in% c('T') & !is.na(pdl1)) %>%
  mutate(pdl1 = case_when(pdl1 %in% c(0,'<1') ~ 'neg', pdl1 %in% c('1-50','>50') ~ 'pos'),
         group7 = paste(treatment, pdl1, sample_time, sep = '.'))
meta = meta %>% dplyr::select(c('sample_id','patient_id','group7'))
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% dplyr::select(-subtype)
subtype = merge(subtype, meta, by = 'sample_id') %>%
  mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1', immune == 'Immune-Enriched, Fibrotic' ~ 'C2', immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'))
# Visualization
prop = subtype %>%
  group_by(group7, immune) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group7) %>%
  mutate(prop = count / sum(count),
         percent = paste0(round(prop * 100, 1), '%'),
         group7 = factor(group7, levels = c('combo.pos.pre','combo.pos.post','combo.neg.pre','combo.neg.post',
                                            'single.pos.pre','single.pos.post','single.neg.pre','single.neg.post'))) %>% ungroup()
ggplot(prop, aes(x = group7, y = prop, fill = immune)) +
  geom_col(position = 'stack', width = 0.7, color = 'white', size = 0.2) +
  geom_label(aes(label = percent), position = position_stack(vjust = 0.5), size = 3.5, color = 'black', fontface = 'bold') +
  scale_fill_manual(values = c('C1'='#ad4445','C2'='#e1b76e','C3'='#9cba59','C4'='#0c6868')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 16, margin = margin(b = 15)),
        plot.subtitle = element_text(hjust = 0.5, color = 'gray40', size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12, margin = margin(r = 10)),
        axis.text.x = element_text(face = 'bold', size = 12, color = 'black'),
        axis.text.y = element_text(size = 11, color = 'black'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(face = 'bold', size = 12),
        legend.text = element_text(size = 11),
        legend.box.background = element_rect(color = 'gray80', fill = 'white'),
        plot.caption = element_text(color = 'gray50', size = 10, hjust = 0))

#### Subtype Sankey plot
# Prepare data
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>%
  dplyr::filter(sample_type %in% c('T') & !is.na(pdl1)) %>%
  mutate(pdl1 = case_when(pdl1 %in% c(0,'<1') ~ 'neg', pdl1 %in% c('1-50','>50') ~ 'pos'))
meta.sample = meta %>% dplyr::select(c('sample_id','patient_id','sample_time'))
meta.patient = meta %>% dplyr::select(c('patient_id','treatment','pdl1')) %>% unique()
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% dplyr::select(-subtype)
subtype = merge(subtype, meta.sample, by = 'sample_id') %>%
  mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1', immune == 'Immune-Enriched, Fibrotic' ~ 'C2', immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'))
# Summarize subtype transitions
pair = subtype %>% pivot_wider(id_cols = patient_id, names_from = sample_time, values_from = immune) %>%
  left_join(meta.patient, by = 'patient_id') %>% na.omit()
# Visualization
pair %>% dplyr::filter(treatment == 'combo' & pdl1 == 'pos') %>%
  group_by(pre, post) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(pre = factor(pre, levels = c('C1','C2','C3','C4'), labels = c('C1b','C2b','C3b','C4b')),
         post = factor(post, levels = c('C1','C2','C3','C4'), labels = c('C1p','C2p','C3p','C4p'))) %>%
ggplot(., aes(axis1 = pre, axis2 = post, y = count)) +
  geom_alluvium(aes(fill = pre), alpha = 0.8, width = 1/6, color = 'white', curve_type = 'cubic', knot.pos = 0.5) +
  geom_stratum(aes(fill = post), alpha = 0.9, width = 1/6, color = 'white') +
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 3.5, color = 'black', fontface = 'bold') +
  geom_text(stat = 'flow', aes(label = ifelse(count > 0, count, '')), size = 3, color = 'white', fontface = 'bold', nudge_x = 0.05) +
  scale_fill_manual(values = c('C1p'='#ae2012','C2p'='#ca6702','C3p'='#0a9396','C4p'='#005f73',
                               'C1b'='#ae2012','C2b'='#ca6702','C3b'='#0a9396','C4b'='#005f73')) +
  scale_x_discrete(limits = c('Pre-treatment', 'Post-treatment'), expand = c(0.05, 0.05)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 16),
        plot.subtitle = element_text(hjust = 0.5, color = 'gray40'),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(face = 'bold'),
        plot.caption = element_text(color = 'gray50', size = 10))

#### Gene expression
# Clinical metadata
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>%
  dplyr::filter(sample_type %in% c('T') & !is.na(pdl1)) %>%
  mutate(patient_id = factor(patient_id),
         sample_time = factor(sample_time, levels = c('pre','post')),
         pdl1 = case_when(pdl1 %in% c(0,'<1') ~ 'neg', pdl1 %in% c('1-50','>50') ~ 'pos'),
         group = paste(treatment, pdl1, sample_time, sep = '.'))
anno = meta %>% dplyr::select(c('sample_id','orr','response','pdl1','group7')) %>%
  merge(., subtype %>% dplyr::select('sample_id','immune'), by = 'sample_id') %>% rename(group = group7) %>%
  mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1', immune == 'Immune-Enriched, Fibrotic' ~ 'C2', immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'),
         group = factor(group, levels = c('combo.pos.pre','combo.pos.post','combo.neg.pre','combo.neg.post','single.pos.pre','single.pos.post','single.neg.pre','single.neg.post')))
# Transcriptome matrix
matrix = read.csv('tpm.csv', header = TRUE, row.names = 1) %>% dplyr::select(meta$sample_id)
matrix = log2(matrix + 1)
gene.immune = matrix %>% dplyr::filter(rownames(.) %in% c('CD8A','CD8B','CD4','PDCD1','CTLA4','LAG3','FGFBP2','FCGR3A','KLRD1','IL2','TNFA','IFNA','IFNG','PRF1','TRDC','TRDV1','TRDV3'))
subtype = read.csv('ccsubtype.csv', header = TRUE)

#### Immune microenvironment features
# Calculate immune pathway scores
innate = read.gmt('signatures-pathway.gmt')
colnames(innate) = c('term','SYMBOL')
innate = lapply(split(innate$SYMBOL, innate$term), function(x) {x = x[x != '']})
sig.innate = gsva(as.matrix(matrix), innate, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
# Calculate immune cell infiltration scores
nature = read.gmt('signatures-cell.gmt')
colnames(nature) = c('term','SYMBOL')
nature = lapply(split(nature$SYMBOL, nature$term), function(x) {x = x[x != '']})
sig.nature = gsva(as.matrix(matrix), nature, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

# Visualize clustering results
anno = anno %>% arrange(group, immune) %>% dplyr::select(group, immune)
df = gene.immune %>% t() %>% scale() %>% t() %>% as.data.frame() %>% dplyr::select(rownames(anno))
gap = cumsum(table(anno$group))[-length(table(anno$group))]
bk = c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(df, annotation_col = anno, gaps_col = gap,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = FALSE,
         legend = FALSE, fontsize_row = 5,
         color = c(colorRampPalette(colors = c('#1e96fc','#f8f9fa'))(length(bk)/2),
                   colorRampPalette(colors = c('#f8f9fa','#f77f00'))(length(bk)/2)),
         annotation_colors = list(
           recist = c('CR'='#e29578', 'PR'='#ffdab9', 'SD'='#83c5be', 'PD'='#006d77'),
           response = c('cpr'='#ed6a5a', 'mpr'='#ffd166', 'non-mpr'='#5ca4a9', 'non-surgery'='#adb5bd'),
           immune = c('C1'='#bc4749', 'C2'='#a7c957', 'C3'='#6a994e', 'C4'='#386641')))


#### Group comparison
# Reshape to long format
compare = gene.immune %>% as.data.frame() %>% rownames_to_column('obj') %>%
  pivot_longer(., -obj, names_to = 'sample_id', values_to = 'expr') %>%
  merge(., anno %>% dplyr::select(c('sample_id','group')), by = 'sample_id')
colnames(compare) = c('sample_id','obj','expr','group')
compare = compare %>% filter(!is.na(group))
compare$obj = factor(compare$obj, levels = c('PDCD1','CTLA4','TIGIT','LAG3','HAVCR2','HLA-E','IFNG','PRF1','TRDC','TRDV1','TRDV2','TRDV3','TLR3','TLR7'))
compare$obj = factor(compare$obj, levels = c('APM','APM-MHCI','APM-MHCII','DC-APM','TC-APM','AIR','IIR','IIRAST','Toll'))
compare$obj = factor(compare$obj, levels = c('B_cells','T_cells','CD8_T_cells','Cytotoxic_cells','GD_1_3_T_cells','NK_cells','Dendritic_cells','Monocytes','Macrophages','Mast_cells','Neutrophils'))

# Differential analysis across groups
res = compare %>% group_by(obj) %>%
  summarise(combo.pos_pval = wilcox.test(expr[group == levels(compare$group)[1]], expr[group == levels(compare$group)[2]], paired = FALSE, exact = FALSE)$p.value,
            combo.neg_pval = wilcox.test(expr[group == levels(compare$group)[3]], expr[group == levels(compare$group)[4]], paired = FALSE, exact = FALSE)$p.value,
            single.pos_pval = wilcox.test(expr[group == levels(compare$group)[5]], expr[group == levels(compare$group)[6]], paired = FALSE, exact = FALSE)$p.value,
            single.neg_pval = wilcox.test(expr[group == levels(compare$group)[7]], expr[group == levels(compare$group)[8]], paired = FALSE, exact = FALSE)$p.value,
            combo.pos_diff = log2(mean(expr[group == levels(compare$group)[2]]) / mean(expr[group == levels(compare$group)[1]])),
            combo.neg_diff = log2(mean(expr[group == levels(compare$group)[4]]) / mean(expr[group == levels(compare$group)[3]])),
            single.pos_diff = log2(mean(expr[group == levels(compare$group)[6]]) / mean(expr[group == levels(compare$group)[5]])),
            single.neg_diff = log2(mean(expr[group == levels(compare$group)[8]]) / mean(expr[group == levels(compare$group)[7]]))
  ) %>%
  mutate(combo.pos_padj = round(p.adjust(combo.pos_pval, method = 'BH'), 10),
         combo.neg_padj = round(p.adjust(combo.neg_pval, method = 'BH'), 10),
         single.pos_padj = round(p.adjust(single.pos_pval, method = 'BH'), 10),
         single.neg_padj = round(p.adjust(single.neg_pval, method = 'BH'), 10))

# Draw bubble plot
plot = res %>% dplyr::select(c('obj','combo.pos_diff':'single.neg_diff','combo.pos_padj':'single.neg_padj'))
plot = plot %>% pivot_longer(cols = -obj, names_to = c('comp', '.value'), names_sep = '_') %>%
  mutate(comp = factor(comp, levels = c('combo.pos', 'combo.neg', 'single.pos', 'single.neg')),
         color = ifelse(diff > 0, -log10(padj), -(-log10(padj))), size = abs(diff)) %>%
  dplyr::filter(!is.infinite(diff) & !is.na(diff) & diff != 'NaN')
summary((plot$size))
summary((plot$color))
ggplot(plot, aes(x = obj, y = comp)) +
  geom_point(aes(size = size, color = color), alpha = 0.9) +
  geom_text(data = subset(plot, padj < 0.1), aes(label = '*'), color = 'black', size = 5, vjust = 0.7) +
  scale_size_continuous(name = 'log2fc', range = c(2,8), breaks = seq(0, 0.4, length.out = 5)) +
  scale_color_gradient2(name = '-log10(adj.P)', low = '#0096c7', mid = 'white', high = '#d00000', midpoint = 0, limits = c(-4, 4)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 6), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = 'grey90', linewidth = 0.15))

# Draw comparison bubble plot
plot = res %>% dplyr::select(c('obj','combo.pos_diff','combo.pos_pval','single.pos_diff','single.pos_pval','combo.neg_diff','combo.neg_pval','single.neg_diff','single.neg_pval')) %>%
  mutate(sig = ifelse(combo.pos_diff > single.pos_diff + .5, 'combo', ifelse(combo.pos_diff < single.pos_diff - 1, 'single', 'none')))
ggplot(plot, aes(x = combo.pos_diff, y = single.pos_diff)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_abline(slope = 1, intercept = 1, linetype = 'dotted', color = 'gray60', alpha = 0.5) +
  geom_abline(slope = 1, intercept = -1, linetype = 'dotted', color = 'gray60', alpha = 0.5) +
  geom_point(aes(color = sig), size = 4, alpha = 0.7) +
  geom_text_repel(data = subset(plot, sig != 'not sig'), aes(label = obj), size = 3, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', segment.size = 0.5, segment.alpha = 0.6, force = 1) +
  # scale_color_manual(name = 'significance', values = c('both' = '#E41A1C', 'combo only' = '#377EB8', 'single only' = '#4DAF4A', 'not sig' = 'gray70')) +
  scale_color_manual(name = 'significance', values = c('combo' = '#E41A1C', 'single' = '#377EB8', 'none' = 'gray60')) +
  scale_x_continuous(limits = c(-1, 5), breaks = seq(-1, 5, 1)) +
  scale_y_continuous(limits = c(-1, 5), breaks = seq(-1, 5, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'grey90', linewidth = 0.15),
        legend.position = 'none') +
  labs(x = 'combo diff', y = 'single diff')

# Draw boxplots
df = meta %>% merge(., gene.immune %>% t() %>% as.data.frame() %>% rownames_to_column('sample_id')) %>%
  merge(., sig.nature %>% t() %>% as.data.frame() %>% rownames_to_column('sample_id'))
df[which(colnames(df) == 'CD4'):ncol(df)] = df[which(colnames(df) == 'CD4'):ncol(df)] %>% scale()
df.combo.pos = df %>% dplyr::filter(treatment == 'combo' & pdl1 == 'pos')
df.combo.neg = df %>% dplyr::filter(treatment == 'combo' & pdl1 == 'neg')
df.single.pos = df %>% dplyr::filter(treatment == 'single' & pdl1 == 'pos')
df.single.neg = df %>% dplyr::filter(treatment == 'single' & pdl1 == 'neg')
# Boxplots for all samples
ggplot(df.combo.pos, aes(sample_time, CD8_T_cells, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = 'grey20') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
ggplot(df.combo.neg, aes(sample_time, CD8_T_cells, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = 'grey20') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
ggplot(df.single.pos, aes(sample_time, CD8_T_cells, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = 'grey20') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
ggplot(df.single.neg, aes(sample_time, CD8_T_cells, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = 'grey20') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
# Boxplots for paired samples
ggplot(df.combo.pos %>% dplyr::filter(patient_id %in% df.combo.pos$patient_id[duplicated(df.combo.pos$patient_id)]),
       aes(sample_time, TRDV1, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')
ggplot(df.combo.neg %>% dplyr::filter(patient_id %in% df.combo.neg$patient_id[duplicated(df.combo.neg$patient_id)]),
       aes(sample_time, TRDV1, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')
ggplot(df.single.pos %>% dplyr::filter(patient_id %in% df.single.pos$patient_id[duplicated(df.single.pos$patient_id)]),
       aes(sample_time, TRDV1, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')
ggplot(df.single.neg %>% dplyr::filter(patient_id %in% df.single.neg$patient_id[duplicated(df.single.neg$patient_id)]),
       aes(sample_time, TRDV1, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')


#### Paired-sample comparison
# Clinical metadata
dup = meta$patient_id[duplicated(meta$patient_id)]
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% dplyr::select('sample_id','immune')
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>%
  mutate(patient_id = factor(patient_id), sample_time = factor(sample_time, levels = c('pre','post')),
         pdl1 = case_when(pdl1 %in% c(0,'<1') ~ 'neg', pdl1 %in% c('1-50','>50') ~ 'pos'),
         group = paste(treatment, pdl1, sample_time, sep = '.')) %>%
  dplyr::filter(patient_id %in% dup) %>%
  dplyr::select(c('sample_id','treatment','orr','response','pdl1','group')) %>%
  merge(., subtype, by = 'sample_id') %>%
  mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1', immune == 'Immune-Enriched, Fibrotic' ~ 'C2', immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'),
         group = factor(group, levels = c('combo.pos.pre','combo.pos.post','combo.neg.pre','combo.neg.post','single.pos.pre','single.pos.post','single.neg.pre','single.neg.post')))
meta.combo.pos = meta %>% dplyr::filter(treatment == 'combo' & pdl1 == 'pos')
meta.combo.neg = meta %>% dplyr::filter(treatment == 'combo' & pdl1 == 'neg')
meta.single.pos = meta %>% dplyr::filter(treatment == 'single' & pdl1 == 'pos')
meta.single.neg = meta %>% dplyr::filter(treatment == 'single' & pdl1 == 'neg')
# Transcriptome matrix
matrix = read.csv('tpm.csv', header = TRUE, row.names = 1) %>% dplyr::select(meta$sample_id)
matrix = log2(matrix + 1)
target = c('CD8A','CD8B','CD4','PDCD1','CTLA4','TIGIT','IFNG','PRF1','TRDC','TRDV1','TRDV3')
gene.combo.pos = matrix %>% dplyr::filter(rownames(.) %in% target) %>% dplyr::select(meta.combo.pos$sample_id)
gene.combo.neg = matrix %>% dplyr::filter(rownames(.) %in% target) %>% dplyr::select(meta.combo.neg$sample_id)
gene.single.pos = matrix %>% dplyr::filter(rownames(.) %in% target) %>% dplyr::select(meta.single.pos$sample_id)
gene.single.neg = matrix %>% dplyr::filter(rownames(.) %in% target) %>% dplyr::select(meta.single.neg$sample_id)

## Immune microenvironment features
# Immune-related pathways
innate = read.gmt('signatures-pathway.gmt')
colnames(innate) = c('term','SYMBOL')
innate = lapply(split(innate$SYMBOL, innate$term), function(x) {x = x[x != '']})
innate.combo.pos = gsva(as.matrix(matrix), innate, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.combo.pos$sample_id)
innate.combo.neg = gsva(as.matrix(matrix), innate, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.combo.neg$sample_id)
innate.single.pos = gsva(as.matrix(matrix), innate, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.single.pos$sample_id)
innate.single.neg = gsva(as.matrix(matrix), innate, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.single.neg$sample_id)
# Immune infiltration cells
nature = read.gmt('signatures-cell.gmt')
colnames(nature) = c('term','SYMBOL')
nature = lapply(split(nature$SYMBOL, nature$term), function(x) {x = x[x != '']})
nature.combo.pos = gsva(as.matrix(matrix), nature, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.combo.pos$sample_id)
nature.combo.neg = gsva(as.matrix(matrix), nature, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.combo.neg$sample_id)
nature.single.pos = gsva(as.matrix(matrix), nature, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.single.pos$sample_id)
nature.single.neg = gsva(as.matrix(matrix), nature, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>% as.data.frame() %>% dplyr::select(meta.single.neg$sample_id)

## Differential analysis across groups
res.combo.pos = nature.combo.pos %>% rownames_to_column('feature') %>%
  pivot_longer(-feature, names_to = 'sample', values_to = 'value') %>%
  separate(sample, into = c('patient_id', 'sample_time'), sep = '\\.') %>%
  pivot_wider(names_from = sample_time, values_from = value) %>%
  mutate(diff = post - pre) %>%
  group_by(feature) %>%
  summarise(combo.pos_pval = wilcox.test(post, pre, paired = TRUE, exact = FALSE)$p.value,
            combo.pos_diff = mean(diff, na.rm = TRUE), .groups = 'drop')
res.combo.neg = nature.combo.neg %>% rownames_to_column('feature') %>%
  pivot_longer(-feature, names_to = 'sample', values_to = 'value') %>%
  separate(sample, into = c('patient_id', 'sample_time'), sep = '\\.') %>%
  pivot_wider(names_from = sample_time, values_from = value) %>%
  mutate(diff = post - pre) %>%
  group_by(feature) %>%
  summarise(combo.neg_pval = wilcox.test(post, pre, paired = TRUE, exact = FALSE)$p.value,
            combo.neg_diff = mean(diff, na.rm = TRUE), .groups = 'drop')
res.single.pos = nature.single.pos %>% rownames_to_column('feature') %>%
  pivot_longer(-feature, names_to = 'sample', values_to = 'value') %>%
  separate(sample, into = c('patient_id', 'sample_time'), sep = '\\.') %>%
  pivot_wider(names_from = sample_time, values_from = value) %>%
  mutate(diff = post - pre) %>%
  group_by(feature) %>%
  summarise(single.pos_pval = wilcox.test(post, pre, paired = TRUE, exact = FALSE)$p.value,
            single.pos_diff = mean(diff, na.rm = TRUE), .groups = 'drop')
res.single.neg = nature.single.neg %>% rownames_to_column('feature') %>%
  pivot_longer(-feature, names_to = 'sample', values_to = 'value') %>%
  separate(sample, into = c('patient_id', 'sample_time'), sep = '\\.') %>%
  pivot_wider(names_from = sample_time, values_from = value) %>%
  mutate(diff = post - pre) %>%
  group_by(feature) %>%
  summarise(single.neg_pval = wilcox.test(post, pre, paired = TRUE, exact = FALSE)$p.value,
            single.neg_diff = mean(diff, na.rm = TRUE), .groups = 'drop')
# Draw comparison bubble plot
plot = merge(res.combo.pos, res.single.pos, by = 'feature') %>% dplyr::select(c('feature','combo.pos_diff','combo.pos_pval','single.pos_diff','single.pos_pval')) %>%
  mutate(sig = ifelse(combo.pos_diff > single.pos_diff + .1, 'combo', 'single'))
plot = merge(res.combo.neg, res.single.neg, by = 'feature') %>% dplyr::select(c('feature','combo.neg_diff','combo.neg_pval','single.neg_diff','single.neg_pval')) %>%
  mutate(sig = ifelse(combo.neg_diff > single.neg_diff + .1, 'combo', 'single'))
ggplot(plot, aes(x = combo.neg_diff, y = single.neg_diff)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray70', alpha = 0.7) +
  geom_abline(slope = 1, intercept = -.5, linetype = 'dotted', color = 'gray60', alpha = 0.5) +
  geom_point(aes(color = sig), size = 4, alpha = 0.7) +
  geom_text_repel(data = subset(plot, sig != 'not sig'), aes(label = feature), size = 3, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', segment.size = 0.5, segment.alpha = 0.6, force = 1) +
  scale_color_manual(name = 'significance', values = c('combo' = '#E41A1C', 'single' = '#377EB8')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'grey90', linewidth = 0.15),
        legend.position = 'none') +
  scale_x_continuous(limits = c(-.2, .4), breaks = seq(-.2, .4, .1)) +
  scale_y_continuous(limits = c(-.2, .4), breaks = seq(-.2, .4, .1)) +
  labs(x = 'combo diff', y = 'single diff')

