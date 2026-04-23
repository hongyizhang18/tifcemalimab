# This script visualizes subtype changes before and after treatment.
# Updated: 2026-01-17

# Load packages
ps = c('tidyverse','openxlsx','clusterProfiler','GSVA','ggsankey','ggplot2','ggalluvial','pheatmap')
for (i in ps) {
  library(i, character.only = TRUE)
}

#### Subtype bar plot
# Prepare data
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>% dplyr::filter(sample_type %in% c('T'))
meta = meta %>% dplyr::select(c('sample_id','patient_id','group','orr','mpr','cpr'))
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% dplyr::select(-subtype)
subtype = merge(subtype, meta, by = 'sample_id') %>%
  mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1',
                            immune == 'Immune-Enriched, Fibrotic' ~ 'C2',
                            immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'))
# Visualize subtype proportions
prop = subtype %>%
  group_by(group, immune) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(prop = count / sum(count),
         percent = paste0(round(prop * 100, 1), '%'),
         group = factor(group, levels = c('combo.T.pre', 'combo.T.post', 'single.T.pre', 'single.T.post'))) %>%
  ungroup()
ggplot(prop, aes(x = group, y = prop, fill = immune)) +
  geom_col(position = 'stack', width = 0.7, color = 'white', size = 0.2) +
  geom_label(aes(label = percent), position = position_stack(vjust = 0.5), size = 3.5, color = 'black', fontface = 'bold') +
  scale_fill_manual(values = c('C1'='#bc4749', 'C2'='#eec170', 'C3'='#a7c957', 'C4'='#0D6E6E')) +
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
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>% dplyr::filter(sample_type %in% c('T'))
meta.sample = meta %>% dplyr::select(sample_id, patient_id, sample_time)
meta.patient = meta %>% dplyr::select(patient_id, treatment, orr, mpr1, cpr1) %>% unique()
subtype = read.csv('ccsubtype.csv', header = TRUERUE) %>% dplyr::select(-subtype)
# Process cluster-based subtypes
subtype = merge(subtype, meta.sample, by = 'sample_id') %>%
  mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1',
                            immune == 'Immune-Enriched, Fibrotic' ~ 'C2',
                            immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'))
pair = subtype %>% 
  pivot_wider(id_cols = patient_id, names_from = sample_time, values_from = immune) %>%
  left_join(meta.patient, by = 'patient_id') %>% 
  na.omit()
# Build Sankey input table
sankey = pair %>% dplyr::filter(treatment == 'combo') %>%
  dplyr::select(patient_id, pre, post) %>%
  make_long(pre, post)
sankey = sankey %>% mutate(node = factor(node, levels = c('C4','C3','C2','C1')),
                           next_node = factor(next_node, levels = c('C4','C3','C2','C1')))
# Visualize with ggsankey
ggplot(sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node)) +
  geom_sankey(flow.alpha = 0.7, node.width = 0.1, node.color = 'white') +
  scale_fill_manual(values = c('C1'='#bc4749', 'C2'='#eec170', 'C3'='#a7c957', 'C4'='#0D6E6E')) +
  scale_x_discrete(labels = c('Pre-treatment', 'Post-treatment')) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')


#### Subtype Sankey plot
# Prepare data
meta = read.xlsx('sample-info.xlsx', sheet = 1) %>% dplyr::filter(sample_type %in% c('T'))
meta.sample = meta %>% dplyr::select(c('sample_id','patient_id','sample_time'))
meta.patient = meta %>% dplyr::select(c('patient_id','treatment','orr','mpr','cpr')) %>% unique()
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% dplyr::select(-subtype)
subtype = merge(subtype, meta.sample, by = 'sample_id') %>%
  mutate(immune = case_when(immune == 'Immune-Enriched, Non-Fibrotic' ~ 'C1',
                            immune == 'Immune-Enriched, Fibrotic' ~ 'C2',
                            immune == 'Fibrotic' ~ 'C3', immune == 'Depleted' ~ 'C4'))
# Visualize Sankey plot
pair = subtype %>% pivot_wider(id_cols = patient_id, names_from = sample_time, values_from = immune) %>%
  left_join(meta.patient, by = 'patient_id') %>% na.omit()
pair = pair %>%
  dplyr::filter(treatment == 'single') %>%
  group_by(pre, post) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(pre = factor(pre, levels = c('C1','C2','C3','C4'), labels = c('C1b','C2b','C3b','C4b')),
         post = factor(post, levels = c('C1','C2','C3','C4'), labels = c('C1p','C2p','C3p','C4p')))
ggplot(pair, aes(axis1 = pre, axis2 = post, y = count)) +
  geom_alluvium(aes(fill = pre), alpha = 0.8, width = 1/8, color = 'white', curve_type = 'cubic', knot.pos = 0.5) +
  geom_stratum(aes(fill = post), alpha = 0.9, width = 1/8, color = 'white') +
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 3.5, color = 'black', fontface = 'bold') +
  geom_text(stat = 'flow', aes(label = ifelse(count > 0, count, '')), size = 3, color = 'white', fontface = 'bold', nudge_x = 0.05) +
  scale_fill_manual(values = c('C1p'='#bc4749','C2p'='#eec170','C3p'='#a7c957','C4p'='#0D6E6E',
                               'C1b'='#bc4749','C2b'='#eec170','C3b'='#a7c957','C4b'='#0D6E6E')) +
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

