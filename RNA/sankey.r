# This script draws a Sankey diagram for patient transitions.
# Updated: 2025-02-10

# Load packages
ps = c('tidyverse','openxlsx','clusterProfiler','GSVA','ggsankey','ggplot2','ggalluvial','networkD3','htmlwidgets')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Data import
# Transcriptome matrix
tpm = read.csv('tpm.csv', header = TRUE, row.names = 1)
tpm = log2(tpm + 1)
gene = tpm %>% dplyr::filter(rownames(.) %in% c('TRDC','TRDV1','TRDV3')) %>% t() %>% as.data.frame()
# Clinical metadata
dup = read.xlsx('sample-info.xlsx', sheet = 1) %>% dplyr::filter(sample_type == 'T') %>% pull(patient_id) %>% .[duplicated(.)]
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>%
  dplyr::filter(sample_type == 'T' & treatment == 'combo' & sample_time == 'pre' & !is.na(pdl1))
gene = gene %>% dplyr::filter(rownames(.) %in% meta$sample_id)
rownames(gene) = str_sub(rownames(gene), 1, 5)
# Integrate data
sig = gsva(as.matrix(tpm), list(gdt1 = c('TRDV1','TRDV3'), gdt2 = c('TRDC','TRDV1','TRDV3')), method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE) %>%
  t() %>% as.data.frame() %>% scale() %>% as.data.frame() %>% dplyr::filter(grepl('pre', rownames(.)))
rownames(sig) = str_sub(rownames(sig), 1, 5)
df = meta %>% merge(., sig %>% rownames_to_column('patient_id'), by = 'patient_id') %>%
  merge(., gene %>% rownames_to_column('patient_id'), by = 'patient_id') %>%
  dplyr::select(c('patient_id','recist','response','pdl1','gdt1','gdt2'))
# 更改数据
factor.group = c('neg.lo','neg.hi','pos.lo','pos.hi')
factor.recist = c('PD','SD','PR','CR')
factor.response = c('non-surgery due to progressive disease',
                    'non-surgery due to patients’ decision',
                    'non-surgery due to surgical infeasibility',
                    'non-mpr','non-surgery due to tumor lesion disappearance','mpr','cpr')
df = df %>% mutate(pdl1_group = ifelse(pdl1 %in% c(0,'<1'), 'neg', 'pos'),
                   gdt2_group = ifelse(gdt2 > median(.$gdt2), 'hi', 'lo'),
                   group = paste(pdl1_group, gdt2_group, sep = '.'))

#### 绘制桑基图1
# 预处理
sankey = df %>% make_long(recist, group, response) %>%
  mutate(node = factor(node, levels = c(factor.group, factor.recist, factor.response)),
         next_node = factor(next_node, levels = c(factor.group, factor.recist, factor.response)))
ggplot(sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node))) +
  geom_sankey() +
  scale_fill_manual(values = c(
    'pos.hi'='#da696c', 'pos.lo'='#edafb8','neg.hi'='#dedbd2', 'neg.lo'='#84a59d',
    'CR' = '#b23a48', 'PR' = '#f4a261', 'SD' = '#83c5be', 'PD' = '#006d77',
    'cpr' = '#d9ed92', 'mpr' = '#b5e48c', 'non-mpr' = '#76c893',
    'non-surgery due to tumor lesion disappearance' = '#99d98c',
    'non-surgery due to surgical infeasibility' = '#52b69a',
    'non-surgery due to patients’ decision' = '#34a0a4',
    'non-surgery due to progressive disease' = '#168aad')) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')

sankey = df %>% group_by(group, recist, response) %>% summarise(count = n(), .groups = 'drop')
# Visualization
ggplot(sankey, aes(axis1 = group, axis2 = recist, axis3 = response, y = count)) +
  geom_alluvium(aes(fill = group), alpha = 0.8, width = 1/6, color = 'white', curve_type = 'cubic', knot.pos = 0.5) +
  geom_stratum(alpha = 0.9, width = 1/6, color = 'white') +
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 3.5, color = 'black', fontface = 'bold') +
  geom_text(stat = 'flow', aes(label = ifelse(count > 0, count, '')), size = 3, color = 'white', fontface = 'bold', nudge_x = 0.05) +
  scale_fill_manual(values = c('pos.hi'='#da696c', 'pos.lo'='#edafb8','neg.hi'='#dedbd2', 'neg.lo'='#84a59d')) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')

