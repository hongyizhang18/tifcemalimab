# merge mutation data with clinical features and derive TMB and ITH metrics.

packages <- c("tidyverse", "openxlsx", "data.table", "maftools")
invisible(lapply(packages, library, character.only = TRUE))

# Read clinical metadata.
sample <- read.csv("WES-metadata.csv", header = TRUE) %>%
  dplyr::filter(Tumor != "") %>%
  dplyr::filter(!grepl("p$", Patient))

clinic <- read.xlsx("patient-info.xlsx", sheet = 1) %>%
  dplyr::filter(patient_id %in% sample$Patient) %>%
  dplyr::select(patient_id, treatment, recist, response) %>%
  mutate(
    treatment = ifelse(treatment == 2, "combo", "single"),
    cpr = case_when(
      response == "CPR" ~ "cpr",
      response == "Non-surgery" ~ NA_character_,
      TRUE ~ "ncpr"
    )
  )
colnames(clinic) <- c("Tumor_Sample_Barcode", "treatment", "recist", "response", "cpr")

# Read mutation and genomic signatures.
load("vep.1000g.clean.rdata")
maf <- read.maf(maf.data, isTCGA = FALSE)

hrd <- read.table("scarhrd.txt", header = TRUE)
colnames(hrd) <- c("Tumor_Sample_Barcode", "HRD", "ALT", "LST", "HRD.sum")

msi <- read.table("msisensor2.txt", header = FALSE)
colnames(msi) <- c("Tumor_Sample_Barcode", "MSI")

# Calculate tumor mutational burden.
tmb <- maftools::tmb(maf, captureSize = 35.8) %>%
  dplyr::select(Tumor_Sample_Barcode, total_perMB, total_perMB_log)
colnames(tmb) <- c("Tumor_Sample_Barcode", "TMB", "TMB_log")
tmb <- tmb %>%
  mutate(TMB_group = ifelse(TMB >= 10, "hi", "lo"))
save(tmb, file = "tmb.rdata")

# Infer intratumor heterogeneity from VAF.
maf@data$t_vaf <- maf@data$t_alt_count / maf@data$t_depth
heter <- inferHeterogeneity(
  maf = maf,
  tsb = maf@data$Tumor_Sample_Barcode,
  vafCol = "t_vaf"
)
print(heter$clusterMeans)

# Extract MATH and mean VAF summaries.
math <- heter[["clusterData"]] %>%
  dplyr::select(Tumor_Sample_Barcode, MATH, MedianAbsoluteDeviation) %>%
  unique() %>%
  dplyr::filter(Tumor_Sample_Barcode %in% unique(maf.data$Tumor_Sample_Barcode))

ith <- heter[["clusterMeans"]] %>%
  dplyr::select(Tumor_Sample_Barcode, meanVaf) %>%
  dplyr::filter(Tumor_Sample_Barcode %in% unique(maf.data$Tumor_Sample_Barcode))

save(math, ith, file = "heterogeneity.rdata")

# Plot the MATH distribution.
math_plot <- math %>%
  arrange(MATH) %>%
  rownames_to_column("id") %>%
  mutate(id = as.numeric(id))

ggplot(math_plot, aes(x = id, y = MATH)) +
  geom_point(color = "#999999", size = 1.5) +
  scale_y_continuous(limits = c(20, 120), breaks = seq(20, 120, 20)) +
  geom_hline(yintercept = 67.65, color = "#FF6633", linetype = "dashed", linewidth = 0.8) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "#EEEEEE", linetype = "dashed"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold")
  )
