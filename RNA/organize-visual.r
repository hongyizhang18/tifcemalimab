# This script organizes BTLA sample matching information.
# Updated: 2025-12-14

# Load packages
ps = c('openxlsx','tidyverse','pheatmap')
for (i in ps) {
  library(i, character.only = TRUE)
}

#### Import sample information
# Transcriptome matrix
rna = read.xlsx('sample-rna.xlsx', sheet = 1)
rna = rna %>% dplyr::select(c('patient_id','sample_id','sample_group'))
table(rna$sample_group)
# Genomic samples
dna = read.xlsx('sample-wes.xlsx', sheet = 1)
dna = dna %>% dplyr::select(c('patient_id','sample_id','sample_group'))
table(dna$sample_group)
# Clinical metadata
clinic = read.xlsx('clinic.xlsx', sheet = 1, detectDates = TRUE)
colnames(clinic) = c('patient_id','patient_name','gender','age','smoker','pathology','t','n','stage','treatment',
                     'pdl1','radiology','recist','orr','dcr','residual','response','yn','stage.reduce',
                     'followup','start.time','PFS','PFS.position','PFS.time','OS','OS.reason','OS.time')
clinic = clinic %>% mutate(across(patient_id:followup, ~ ifelse(.x %in% c('/','?','？','unknown','NE','未入组'), NA, .x))) %>%
  mutate(patient_id = gsub('LA-', '', patient_id), radiology = round(as.numeric(radiology), 1))

#### Reshape sample information
# RNA samples
rna.patient = rna %>% arrange(patient_id, sample_group) %>%
  pivot_wider(id_cols = patient_id, names_from = sample_group, values_from = sample_id, names_glue = '{.value}_{sample_group}') %>%
  dplyr::select(c('patient_id','sample_id_T.pre','sample_id_T.post','sample_id_LN.pre','sample_id_LN.post'))
colnames(rna.patient) = c('patient_id','RNA_T.pre','RNA_T.post','RNA_LN.pre','RNA_LN.post')
# DNA samples
dna.patient = dna %>% arrange(patient_id, sample_group) %>%
  pivot_wider(id_cols = patient_id, names_from = sample_group, values_from = sample_id, names_glue = '{.value}_{sample_group}') %>%
  dplyr::select(c('patient_id','sample_id_T.pre','sample_id_N.post'))
colnames(dna.patient) = c('patient_id','DNA_T.pre','DNA_N.post')
# Merge sample tables
sample = merge(rna.patient, dna.patient, by = 'patient_id', all = TRUERUE)
sample = merge(sample, clinic, by = 'patient_id', all.x = TRUERUE)
# Save output data
filter = clinic %>% filter(!patient_name %in% sample$patient_name)
write.xlsx(sample, file = 'patient-info.xlsx')

#### Visualize heatmap(临床组)
# Prepare data
sample = read.xlsx('patient-info.xlsx', sheet = 1, detectDates = TRUE)
matrix = sample %>% dplyr::select(c(1:3,6:7)) %>% column_to_rownames('patient_id') %>% t() %>% as.data.frame()
matrix = as.data.frame(ifelse(is.na(matrix), 0, 1))
anno = sample %>% dplyr::select(-(contains('RNA')|contains('DNA')|contains('S.'))) %>%
  select(-c('patient_name','t','n','radiology','orr','dcr','residual','yn','followup','stage.reduce','start.time')) %>%
  mutate(gender = ifelse(gender == '男', 'male', 'female'),
         age = ifelse(age >= 65, '>=65', '<65'),
         pathology = case_when(grepl('非小', pathology) ~ 'Non-squamous', grepl('腺癌', pathology) ~ 'Non-squamous', grepl('鳞癌', pathology) ~ 'Squamous'),
         treatment = case_when(treatment == 1 ~ 'single', treatment == 2 ~ 'combo'),
         pdl1 = case_when(pdl1 == '0' ~ '<1', pdl1 == '>50%' ~ '>50', TRUE ~ pdl1),
         PFS = case_when(PFS == 1 ~ 'pfs', PFS == 0 ~ 'none', is.na(PFS) ~ 'unknown'),
         OS = case_when(OS == 1 ~ 'os', OS == 0 ~ 'none', is.na(OS) ~ 'unknown')) %>%
  mutate(across(where(is.character), ~ ifelse(is.na(.), 'unknown', .))) %>%
  mutate(gender = factor(gender, levels = c('male','female')),
         age = factor(age, levels = c('<65','>=65')),
         smoker = factor(smoker, levels = c('current','former','never')),
         stage = factor(stage, levels = c('IIIA','IIIB','IIIC')),
         pathology = factor(pathology, levels = c('Squamous','Non-squamous')),
         pdl1 = factor(pdl1, levels = c('<1','1-50','>50','unknown')),
         treatment = factor(treatment, levels = c('combo','single')),
         recist = factor(recist, levels = c('CR','PR','SD','PD')),
         response = factor(response, levels = c('CPR','MPR','Non-MPR','Non-surgery')),
         PFS = factor(PFS, levels = c('pfs','none','unknown')),
         OS = factor(OS, levels = c('os','none','unknown')))
# Reorder annotations
anno = anno %>% arrange(treatment, recist, response) %>% column_to_rownames(var = 'patient_id') %>% dplyr::select(c('treatment','age','gender','smoker','pathology','stage','pdl1','recist','response','PFS','OS'))
matrix = matrix %>% dplyr::select(rownames(anno))
test = matrix %>% t() %>% as.data.frame() %>% dplyr::filter(RNA_T.pre == '0' & RNA_T.post == '0' & DNA_T.pre == '0')
# Draw heatmap
pheatmap(
  matrix, annotation_col = anno,
  cluster_rows = FALSEALSE, cluster_cols = FALSEALSE,
  show_rownames = TRUERUE, show_colnames = FALSEALSE,
  color = colorRampPalette(colors = c('#FFFFFF','#748cab'))(100),
  annotation_colors = list(
    age = c('>=65'='#bcb7ca', '<65'='#f0efeb'),
    gender = c('female'='#f4acb7', 'male'='#f3dad8'),
    smoker = c('current'='#e8871e', 'former'='#ffbf69', 'never'='#f8dda4'),
    pathology = c('Squamous'='#cb997e', 'Non-squamous'='#ede0d4'),
    stage = c('IIIA'='#d8f3dc', 'IIIB'='#76c893', 'IIIC'='#34a0a4'),
    pdl1 = c('<1'='#a8dadc', '1-50'='#5fa8d3', '>50'='#1e6091', 'unknown'='#ced4da'),
    treatment = c('combo'='#ef476f', 'single'='#06d6a0'),
    recist = c('CR'='#e56b6f', 'PR'='#eaac8b', 'SD'='#b56576', 'PD'='#6d597a'),
    response = c('CPR'='#ed6a5a', 'MPR'='#ffd166', 'Non-MPR'='#5ca4a9', 'Non-surgery'='#adb5bd'),
    PFS = c('pfs'='#757bc8', 'none'='#edf2fb', 'unknown'='#ced4da'),
    OS = c('os'='#757bc8', 'none'='#edf2fb', 'unknown'='#ced4da')),
  fontsize = 4,
  fontsize_col = 5,
  fontsize_row = 8,
  border = TRUE)

#### Visualize heatmap(基因组)
# Clinical metadata
sample = read.xlsx('patient-info.xlsx', sheet = 1, detectDates = TRUE)
matrix = sample %>% dplyr::select(c(1:7)) %>% column_to_rownames('patient_id') %>% t() %>% as.data.frame()
matrix = as.data.frame(ifelse(is.na(matrix), 0, 1))
# Multi-omics data
subtype = read.csv('ccsubtype.csv', header = TRUE) %>% dplyr::select(c('sample_id','immune'))
subtype = subtype %>% dplyr::filter(grepl('pre', sample_id)) %>% mutate(sample_id = str_sub(sample_id, 1, 5))
colnames(subtype) = c('patient_id','immune')
hrd = read.table('scarhrd.txt', header = TRUE) %>% dplyr::select(c('X','HRD.sum'))
colnames(hrd) = c('patient_id','hrd')
msi = read.table('msisensor2.txt', header = FALSE) %>% dplyr::filter(!grepl('p$', V1))
colnames(msi) = c('patient_id','msi')
# Integrate data
anno = sample %>% dplyr::select(-(contains('RNA')|contains('DNA')|contains('S.'))) %>%
  dplyr::select(-c('patient_name','t','n','radiology','orr','dcr','residual','yn','followup','stage.reduce','start.time','PFS','OS')) %>%
  mutate(gender = ifelse(gender == '男', 'male', 'female'),
         age = ifelse(age >= 65, '>=65', '<65'),
         pathology = case_when(grepl('非小', pathology) ~ 'NSCLC',
                               grepl('腺癌', pathology) ~ 'LUAD',
                               grepl('鳞癌', pathology) ~ 'LUSC'),
         treatment = case_when(treatment == 1 ~ 'single', treatment == 2 ~ 'combo'),
         pdl1 = gsub('%', '', pdl1),
         recist = factor(recist, levels = c('CR','PR','SD','PD')),
         response = factor(response, levels = c('CPR','MPR','Non-MPR','Non-surgery')))
anno = anno %>% mutate(across(where(is.character), ~ ifelse(is.na(.), 'unknown', .)))
anno = anno %>% merge(., subtype, by = 'patient_id', all.x = TRUE) %>% merge(., hrd, by = 'patient_id', all.x = TRUE) %>% merge(., msi, by = 'patient_id', all.x = TRUE)
# Reorder annotations
anno = anno %>% arrange(treatment, response) %>% column_to_rownames(var = 'patient_id')
matrix = matrix %>% dplyr::select(rownames(anno))
# Draw heatmap
pheatmap(
  matrix, annotation_col = anno,
  cluster_rows = FALSEALSE, cluster_cols = FALSEALSE,
  show_rownames = TRUERUE, show_colnames = FALSEALSE,
  color = colorRampPalette(colors = c('#FFFFFF','#748cab'))(100),
  annotation_colors = list(
    age = c('>=65'='#dedbd2', '<65'='#f7e1d7'),
    gender = c('male'='#b0c4b1', 'female'='#edafb8'),
    smoker = c('current'='#bf616a', 'former'='#d08770', 'never'='#a3be8c'),
    stage = c('IIIA'='#76c893', 'IIIB'='#34a0a4', 'IIIC'='#1a759f'),
    pathology = c('LUAD'='#a3cef1', 'LUSC'='#cad2c5', 'NSCLC'='#e0e1dd'),
    pdl1 = c('0'='#a9d6e5', '<1'='#89c2d9', '1-50'='#61a5c2', '>50'='#2c7da0', 'unknown'='#b3b6b7'),
    treatment = c('combo'='#ef476f', 'single'='#06d6a0'),
    recist = c('CR'='#e56b6f', 'PR'='#eaac8b', 'SD'='#b56576', 'PD'='#6d597a'),
    response = c('CPR'='#ed6a5a', 'MPR'='#ffd166', 'Non-MPR'='#5ca4a9', 'Non-surgery'='#adb5bd'),
    immune = c('Immune-Enriched, Non-Fibrotic'='#bc4749', 'Immune-Enriched, Fibrotic'='#a7c957', 'Fibrotic'='#6a994e', 'Depleted'='#386641')),
  fontsize = 8,
  fontsize_col = 5,
  fontsize_row = 8,
  border = TRUE)

