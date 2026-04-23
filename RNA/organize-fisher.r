# This script evaluates the balance of BTLA patient clinical features.
# Updated: 2025-11-14

# Load packages
ps = c('openxlsx','tidyverse')
for (i in ps) {
  library(i, character.only = TRUE)
}

# Set shared input and output directories.
data_dir <- '.'
output_dir <- '.'

#### Data preprocessing
# Recode variables
df = read.xlsx('patient-info.xlsx', sheet = 1, detectDates = TRUE)
df = df %>% dplyr::select(c('patient_id','treatment','age','gender','smoker','pathology','t','n','stage','pdl1','recist','response')) %>%
  mutate(treatment = case_when(treatment == 1 ~ 'single', treatment == 2 ~ 'combo'),
         age = ifelse(age >= 65, '>=65', '<65'),
         gender = ifelse(gender == '男', 'male', 'female'),
         pathology = case_when(grepl('非小', pathology) ~ 'NSCLC', grepl('腺癌', pathology) ~ 'LUAD', grepl('鳞癌', pathology) ~ 'LUSC'),
         pdl1 = case_when(pdl1 == 0 ~ '<1', pdl1 == '>50%' ~ '>50', is.na(pdl1) ~ 'unknown', TRUE ~ pdl1))
# Summarize distributions
table(df$treatment)
table(df$age, df$treatment)
table(df$gender, df$treatment)
table(df$smoker, df$treatment)
table(df$pathology, df$treatment)
table(df$t, df$treatment)
table(df$n, df$treatment)
table(df$stage, df$treatment)
table(df$pdl1, df$treatment)
table(df$recist, df$treatment)
table(df$response, df$treatment)

#### Contingency table tests
tab = table(df$treatment, df$age)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$gender)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$smoker)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$pathology)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$t)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$n)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$stage)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$pdl1)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$recist)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)
tab = table(df$treatment, df$response)
if (any(chisq.test(tab)$expected < 5) || sum(tab) < 40) fisher.test(tab) else chisq.test(tab)

