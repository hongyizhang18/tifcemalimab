# This script performs gene set scoring for transcriptomic data.
# Updated: 2025-09-27

# Load packages
ps = c('tidyverse','openxlsx','clusterProfiler','GSVA','ConsensusClusterPlus','pheatmap','ggplot2','ggrepel')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Data import
# Genomic samples
load(file.path(data_dir, 'tmb.rdata'))
colnames(tmb) = c('patient_id','TMB','TMB_log','TMB_group')
# Clinical metadata
clinic = read.xlsx('sample-info.xlsx', sheet = 1) %>% dplyr::filter(sample_type %in% c('T')) %>% mutate(sample_time = factor(sample_time, levels = c('pre','post')))
anno = clinic %>% dplyr::select(c('sample_id','orr','response','group'))
anno$group = factor(anno$group, levels = c('combo.r.pre','combo.r.post','combo.nr.pre','combo.nr.post','single.r.pre','single.r.post','single.nr.pre','single.nr.post'))
anno = anno %>% mutate(patient_id = str_sub(sample_id, 1, 5)) %>% merge(tmb, ., by = 'patient_id') %>% mutate(group = paste(group, TMB_group, sep = '.'))
# Transcriptome matrix
matrix = read.csv('tpm.csv', header = TRUE, row.names = 1) %>% dplyr::select(anno$sample_id)
matrix = log2(matrix + 1)
gene.immune = matrix %>% dplyr::filter(rownames(.) %in% c(
  'PDCD1','CTLA4','TIGIT','LAG3','HAVCR2','HLA-E',
  'IFNG','PRF1','TRDC','TRDV1','TRDV2','TRDV3','TLR3','TLR7'))
# Immune infiltration subtype
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% dplyr::filter(sample_id %in% anno$sample_id)


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

#### Group comparison
# Reshape to long format
compare = gene.display %>% as.data.frame() %>% rownames_to_column('obj') %>%
  pivot_longer(., -obj, names_to = 'sample_id', values_to = 'expr') %>%
  merge(., anno %>% dplyr::select(c('sample_id','group')), by = 'sample_id')
colnames(compare) = c('sample_id','obj','expr','group')
compare = compare %>% filter(!is.na(group))
compare$group = factor(compare$group, levels = c(
  'combo.r.pre.hi','combo.r.post.hi','combo.r.pre.lo','combo.r.post.lo',
  'combo.nr.pre.hi','combo.nr.post.hi','combo.nr.pre.lo','combo.nr.post.lo',
  'single.r.pre.hi','single.r.post.hi','single.r.pre.lo','single.r.post.lo',
  'single.nr.pre.hi','single.nr.post.hi','single.nr.pre.lo','single.nr.post.lo'))
compare$obj = factor(compare$obj, levels = c('PDCD1','CTLA4','TIGIT','LAG3','HAVCR2','HLA-E','IFNG','PRF1','TRDC','TRDV1','TRDV2','TRDV3','TLR3','TLR7'))
compare$obj = factor(compare$obj, levels = c('APM','APM-MHCI','APM-MHCII','DC-APM','TC-APM','AIR','IIR','IIRAST','Toll'))
compare$obj = factor(compare$obj, levels = c('B_cells','T_cells','CD8_T_cells','Cytotoxic_cells','GD_1_3_T_cells','NK_cells','Dendritic_cells','Monocytes','Macrophages','Mast_cells','Neutrophils'))

# Differential analysis across groups(差值)
res = compare %>% group_by(obj) %>%
  summarise(
    combo.r.hi_pval = wilcox.test(expr[group == levels(compare$group)[1]], expr[group == levels(compare$group)[2]], paired = FALSE, exact = FALSE)$p.value,
    combo.r.lo_pval = wilcox.test(expr[group == levels(compare$group)[3]], expr[group == levels(compare$group)[4]], paired = FALSE, exact = FALSE)$p.value,
    combo.nr.hi_pval = wilcox.test(expr[group == levels(compare$group)[5]], expr[group == levels(compare$group)[6]], paired = FALSE, exact = FALSE)$p.value,
    combo.nr.lo_pval = wilcox.test(expr[group == levels(compare$group)[7]], expr[group == levels(compare$group)[8]], paired = FALSE, exact = FALSE)$p.value,
    single.r.hi_pval = wilcox.test(expr[group == levels(compare$group)[9]], expr[group == levels(compare$group)[10]], paired = FALSE, exact = FALSE)$p.value,
    single.r.lo_pval = wilcox.test(expr[group == levels(compare$group)[11]], expr[group == levels(compare$group)[12]], paired = FALSE, exact = FALSE)$p.value,
    single.nr.hi_pval = wilcox.test(expr[group == levels(compare$group)[13]], expr[group == levels(compare$group)[14]], paired = FALSE, exact = FALSE)$p.value,
    single.nr.lo_pval = wilcox.test(expr[group == levels(compare$group)[15]], expr[group == levels(compare$group)[16]], paired = FALSE, exact = FALSE)$p.value,
    combo.r.hi_diff = (mean(expr[group == levels(compare$group)[2]]) - mean(expr[group == levels(compare$group)[1]])),
    combo.r.lo_diff = (mean(expr[group == levels(compare$group)[4]]) - mean(expr[group == levels(compare$group)[3]])),
    combo.nr.hi_diff = (mean(expr[group == levels(compare$group)[6]]) - mean(expr[group == levels(compare$group)[5]])),
    combo.nr.lo_diff = (mean(expr[group == levels(compare$group)[8]]) - mean(expr[group == levels(compare$group)[7]])),
    single.r.hi_diff = (mean(expr[group == levels(compare$group)[10]]) - mean(expr[group == levels(compare$group)[9]])),
    single.r.lo_diff = (mean(expr[group == levels(compare$group)[12]]) - mean(expr[group == levels(compare$group)[11]])),
    single.nr.hi_diff = (mean(expr[group == levels(compare$group)[14]]) - mean(expr[group == levels(compare$group)[13]])),
    single.nr.lo_diff = (mean(expr[group == levels(compare$group)[16]]) - mean(expr[group == levels(compare$group)[15]]))
  ) %>%
  mutate(combo.r.hi_padj = round(p.adjust(combo.r.hi_pval, method = 'BH'), 10),
         combo.r.lo_padj = round(p.adjust(combo.r.lo_pval, method = 'BH'), 10),
         combo.nr.hi_padj = round(p.adjust(combo.nr.hi_pval, method = 'BH'), 10),
         combo.nr.lo_padj = round(p.adjust(combo.nr.lo_pval, method = 'BH'), 10),
         single.r.hi_padj = round(p.adjust(single.r.hi_pval, method = 'BH'), 10),
         single.r.lo_padj = round(p.adjust(single.r.lo_pval, method = 'BH'), 10),
         single.nr.hi_padj = round(p.adjust(single.nr.hi_pval, method = 'BH'), 10),
         single.nr.lo_padj = round(p.adjust(single.nr.lo_pval, method = 'BH'), 10))
# Differential analysis across groups(log2fc)
res = compare %>% group_by(obj) %>%
  summarise(
    combo.r.hi_pval = wilcox.test(expr[group == levels(compare$group)[1]], expr[group == levels(compare$group)[2]], paired = FALSE, exact = FALSE)$p.value,
    combo.r.lo_pval = wilcox.test(expr[group == levels(compare$group)[3]], expr[group == levels(compare$group)[4]], paired = FALSE, exact = FALSE)$p.value,
    combo.nr.hi_pval = wilcox.test(expr[group == levels(compare$group)[5]], expr[group == levels(compare$group)[6]], paired = FALSE, exact = FALSE)$p.value,
    combo.nr.lo_pval = wilcox.test(expr[group == levels(compare$group)[7]], expr[group == levels(compare$group)[8]], paired = FALSE, exact = FALSE)$p.value,
    single.r.hi_pval = wilcox.test(expr[group == levels(compare$group)[9]], expr[group == levels(compare$group)[10]], paired = FALSE, exact = FALSE)$p.value,
    single.r.lo_pval = wilcox.test(expr[group == levels(compare$group)[11]], expr[group == levels(compare$group)[12]], paired = FALSE, exact = FALSE)$p.value,
    single.nr.hi_pval = wilcox.test(expr[group == levels(compare$group)[13]], expr[group == levels(compare$group)[14]], paired = FALSE, exact = FALSE)$p.value,
    single.nr.lo_pval = wilcox.test(expr[group == levels(compare$group)[15]], expr[group == levels(compare$group)[16]], paired = FALSE, exact = FALSE)$p.value,
    combo.r.hi_diff = log2(mean(expr[group == levels(compare$group)[2]])) / (mean(expr[group == levels(compare$group)[1]])),
    combo.r.lo_diff = log2(mean(expr[group == levels(compare$group)[4]])) / (mean(expr[group == levels(compare$group)[3]])),
    combo.nr.hi_diff = log2(mean(expr[group == levels(compare$group)[6]])) / (mean(expr[group == levels(compare$group)[5]])),
    combo.nr.lo_diff = log2(mean(expr[group == levels(compare$group)[8]])) / (mean(expr[group == levels(compare$group)[7]])),
    single.r.hi_diff = log2(mean(expr[group == levels(compare$group)[10]])) / (mean(expr[group == levels(compare$group)[9]])),
    single.r.lo_diff = log2(mean(expr[group == levels(compare$group)[12]])) / (mean(expr[group == levels(compare$group)[11]])),
    single.nr.hi_diff = log2(mean(expr[group == levels(compare$group)[14]])) / (mean(expr[group == levels(compare$group)[13]])),
    single.nr.lo_diff = log2(mean(expr[group == levels(compare$group)[16]])) / (mean(expr[group == levels(compare$group)[15]]))
  ) %>%
  mutate(combo.r.hi_padj = round(p.adjust(combo.r.hi_pval, method = 'BH'), 10),
         combo.r.lo_padj = round(p.adjust(combo.r.lo_pval, method = 'BH'), 10),
         combo.nr.hi_padj = round(p.adjust(combo.nr.hi_pval, method = 'BH'), 10),
         combo.nr.lo_padj = round(p.adjust(combo.nr.lo_pval, method = 'BH'), 10),
         single.r.hi_padj = round(p.adjust(single.r.hi_pval, method = 'BH'), 10),
         single.r.lo_padj = round(p.adjust(single.r.lo_pval, method = 'BH'), 10),
         single.nr.hi_padj = round(p.adjust(single.nr.hi_pval, method = 'BH'), 10),
         single.nr.lo_padj = round(p.adjust(single.nr.lo_pval, method = 'BH'), 10))

#### Draw bubble plot
plot = res %>% dplyr::select(c('obj','combo.r.hi_diff':'single.nr.lo_diff','combo.r.hi_padj':'single.nr.lo_padj'))
plot = plot %>% pivot_longer(cols = -obj, names_to = c('comp', '.value'), names_sep = '_') %>%
  mutate(comp = factor(comp, levels = c('combo.r.hi','combo.r.lo','combo.nr.hi','combo.nr.lo','single.r.hi','single.r.lo','single.nr.hi','single.nr.lo')),
         color = ifelse(diff > 0, -log10(padj), -(-log10(padj))), size = abs(diff)) %>%
  dplyr::filter(!is.infinite(diff) & !is.na(diff) & diff != 'NaN')
summary((plot$size))
summary((plot$color))
ggplot(plot, aes(x = obj, y = comp)) +
  geom_point(aes(size = size, color = color), alpha = 0.9) +
  geom_text(data = subset(plot, padj < 0.1), aes(label = '*'), color = 'black', size = 5, vjust = 0.7) +
  scale_size_continuous(name = 'log2fc', range = c(3,8), breaks = seq(0, 0.5, length.out = 6)) +
  scale_color_gradient2(name = '-log10(adj.P)', low = '#0096c7', mid = 'white', high = '#d00000', midpoint = 0, limits = c(-3, 3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 9), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = 'grey90', linewidth = 0.15))

#### Draw boxplots
df = clinic %>% merge(., tmb %>% dplyr::select(c('patient_id','TMB_group')), by = 'patient_id') %>%
  merge(., gene.immune %>% t() %>% as.data.frame() %>% rownames_to_column('sample_id')) %>%
  merge(., sig.nature %>% t() %>% as.data.frame() %>% rownames_to_column('sample_id'))
df[which(colnames(df) == 'CTLA4'):ncol(df)] = df[which(colnames(df) == 'CTLA4'):ncol(df)] %>% scale()
df.combo.hi = df %>% dplyr::filter(treatment == 'combo' & TMB_group == 'hi')
df.combo.lo = df %>% dplyr::filter(treatment == 'combo' & TMB_group == 'lo')
df.single.hi = df %>% dplyr::filter(treatment == 'single' & TMB_group == 'hi')
df.single.lo = df %>% dplyr::filter(treatment == 'single' & TMB_group == 'lo')
# Boxplots for all samples
ggplot(df.combo.hi, aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6, color = 'grey20') +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')
ggplot(df.combo.lo, aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6, color = 'grey20') +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')
ggplot(df.single.hi, aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.4, color = 'grey20') +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')
ggplot(df.single.lo, aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6, color = 'grey20') +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 4) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')
# Boxplots for paired samples
ggplot(df.combo.hi %>% dplyr::filter(patient_id %in% df.combo.hi$patient_id[duplicated(df.combo.hi$patient_id)]),
       aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')
ggplot(df.combo.lo %>% dplyr::filter(patient_id %in% df.combo.lo$patient_id[duplicated(df.combo.lo$patient_id)]),
       aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')
ggplot(df.single.hi %>% dplyr::filter(patient_id %in% df.single.hi$patient_id[duplicated(df.single.hi$patient_id)]),
       aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')
ggplot(df.single.lo %>% dplyr::filter(patient_id %in% df.single.lo$patient_id[duplicated(df.single.lo$patient_id)]),
       aes(sample_time, TRDC, fill = sample_time)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = sample_time), size = 2, alpha = 0.7, color = 'grey20') +
  geom_line(aes(group = patient_id), lwd = 0.8, alpha = 0.6, color = '#69b3a2') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', comparisons = list(c('pre','post')), paired = TRUE) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_fill_manual(values = c('pre' = '#e74c3c', 'post' = '#3498db')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'), panel.grid = element_blank(), legend.position = 'none')

