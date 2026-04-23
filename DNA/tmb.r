# evaluate associations between genomic metrics and clinical response.

packages <- c("tidyverse", "openxlsx", "data.table", "maftools", "ggplot2", "ggpubr")
invisible(lapply(packages, library, character.only = TRUE))

# Read processed genomic metrics.
load("tmb.rdata")
load("heterogeneity.rdata")
load("vep.1000g.clean.rdata")

maf.data <- maf.data %>%
  dplyr::filter(!Tumor_Sample_Barcode %in% c("BT007", "BT016", "BT017", "BT047", "BT052", "BT115p", "BT122p"))
maf <- read.maf(maf.data, isTCGA = FALSE)

hrd <- read.table("scarhrd.txt", header = TRUE)
colnames(hrd) <- c("Tumor_Sample_Barcode", "HRD", "ALT", "LST", "HRD.sum")

msi <- read.table("msisensor2.txt", header = FALSE)
colnames(msi) <- c("Tumor_Sample_Barcode", "MSI")

# Read clinical metadata.
meta <- read.xlsx("clinic.xlsx", sheet = 1) %>%
  dplyr::select("患者编号", "PD", "治疗", "影像学缓解比例（-%）治疗结束", "ORR", "剩余肿瘤", "石蜡报告")
colnames(meta) <- c("Tumor_Sample_Barcode", "PDL1", "treatment", "radiology", "orr", "residual", "response")

# Harmonize clinical labels.
meta <- meta %>%
  mutate(
    Tumor_Sample_Barcode = gsub("LA-", "", Tumor_Sample_Barcode),
    treatment = ifelse(treatment == 1, "single", "combo"),
    orr = case_when(
      orr == "ORR" ~ "r",
      orr == "Non-ORR" ~ "nr",
      orr == "NE" ~ "ne"
    ),
    cpr = case_when(
      response == "CPR" ~ "cpr",
      response %in% c("MPR", "Non-MPR") ~ "ncpr",
      response == "Non-surgery" ~ "ns"
    ),
    mpr = case_when(
      response %in% c("CPR", "MPR") ~ "mpr",
      response == "Non-MPR" ~ "nmpr",
      response == "Non-surgery" ~ "ns"
    ),
    group1 = paste(treatment, orr, sep = "_"),
    group2 = paste(treatment, cpr, sep = "_"),
    group3 = paste(treatment, mpr, sep = "_")
  )

# Combine genomic metrics for correlation analysis.
df <- meta %>%
  merge(tmb, by = "Tumor_Sample_Barcode", all.x = TRUE) %>%
  merge(math, by = "Tumor_Sample_Barcode", all.x = TRUE) %>%
  merge(hrd %>% dplyr::select(Tumor_Sample_Barcode, HRD.sum), by = "Tumor_Sample_Barcode", all.x = TRUE) %>%
  merge(msi, by = "Tumor_Sample_Barcode", all.x = TRUE) %>%
  dplyr::filter(radiology != "NE") %>%
  mutate(
    radiology = as.numeric(radiology),
    residual = as.numeric(na_if(residual, "/"))
  )

ggplot(df %>% dplyr::filter(treatment == "combo"), aes(x = radiology, y = MATH)) +
  geom_point(size = 2, alpha = 0.7, color = "#0081a7") +
  geom_smooth(method = "lm", se = TRUE, color = "#c1121f", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = 0.1, label.y.npc = 0.9, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()

ggplot(df %>% dplyr::filter(treatment == "single"), aes(x = radiology, y = MATH)) +
  geom_point(size = 2, alpha = 0.7, color = "#0081a7") +
  geom_smooth(method = "lm", se = TRUE, color = "#c1121f", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = 0.1, label.y.npc = 0.9, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()

# Test the association between TMB group and ORR.
df <- merge(tmb, meta, by = "Tumor_Sample_Barcode") %>%
  dplyr::filter(radiology != "NE") %>%
  mutate(
    TMB = as.numeric(TMB),
    TMB_log = as.numeric(TMB_log),
    radiology = as.numeric(radiology),
    residual = as.numeric(residual)
  )

df.combo <- df %>%
  dplyr::filter(treatment == "combo")
df.single <- df %>%
  dplyr::filter(treatment == "single")

fisher.test(table(df.combo$TMB_group, df.combo$orr))
fisher.test(table(df.single$TMB_group, df.single$orr))

# Plot ORR proportions by TMB group.
df.prop <- df %>%
  mutate(group = paste(treatment, TMB_group, sep = "_")) %>%
  group_by(group, orr) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(proportion = count / sum(count))

ggplot(df.prop, aes(x = group, y = proportion, fill = orr)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", proportion * 100)),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white",
    fontface = "bold"
  )

# Compare TMB distributions across response groups.
df <- merge(tmb, meta, by = "Tumor_Sample_Barcode") %>%
  dplyr::filter(radiology != "NE") %>%
  mutate(
    group1 = factor(group1, levels = c("combo_r", "combo_nr", "single_r", "single_nr")),
    TMB = as.numeric(TMB),
    TMB_log = as.numeric(TMB_log),
    radiology = as.numeric(radiology),
    residual = as.numeric(residual)
  )

ggplot(df, aes(group1, TMB, fill = group1)) +
  geom_boxplot(width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = c("combo_r" = "#e56b6f", "combo_nr" = "#eaac8b", "single_r" = "#a1c181", "single_nr" = "#619b8a")) +
  stat_compare_means(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p",
    comparisons = combn(levels(as.factor(df$group1)), 2, simplify = FALSE)
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect()
  )
