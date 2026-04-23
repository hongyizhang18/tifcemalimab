# test gene-level mutation associations with pathological response.

packages <- c("tidyverse", "openxlsx", "data.table", "maftools")
invisible(lapply(packages, library, character.only = TRUE))

# Read driver gene annotations.
gene.lusc <- read.xlsx("gene-multiomics.xlsx", sheet = 6, startRow = 2) %>%
  pull(Gene.Symbol) %>%
  unique()

gene.driver <- read.xlsx("gene-multiomics.xlsx", sheet = 5, startRow = 2) %>%
  .[-1, ] %>%
  dplyr::filter(Cancer %in% c("LUAD", "LUSC")) %>%
  pull(Gene) %>%
  unique()

# Read mutation data.
load("vep.1000g.clean.rdata")
patient <- unique(maf.data$Tumor_Sample_Barcode)
length(unique(maf.data$Tumor_Sample_Barcode))

# Read clinical metadata.
sample <- read.csv("WES-metadata.csv", header = TRUE) %>%
  dplyr::filter(Tumor != "") %>%
  dplyr::filter(!grepl("p$", Patient))

clinic <- read.xlsx("patient-info.xlsx", sheet = 1) %>%
  dplyr::filter(patient_id %in% sample$Patient) %>%
  dplyr::select(patient_id, treatment, recist, response) %>%
  mutate(
    treatment = ifelse(treatment == 2, "combo", "single"),
    recist = case_when(
      recist %in% c("CR", "PR") ~ "r",
      recist %in% c("SD", "PD") ~ "nr"
    ),
    cpr = case_when(
      response == "CPR" ~ "cpr",
      response == "Non-surgery" ~ NA_character_,
      TRUE ~ "ncpr"
    )
  )

# Build a binary mutation matrix by patient and gene.
matrix <- data.frame(
  matrix(
    NA,
    nrow = length(patient),
    ncol = length(gene.lusc),
    dimnames = list(patient, gene.lusc)
  )
)

for (i in patient) {
  for (j in colnames(matrix)) {
    gene <- maf.data %>%
      dplyr::filter(Hugo_Symbol == j & Tumor_Sample_Barcode == i)
    matrix[i, j] <- if (nrow(gene) > 0) "mut" else "wild"
  }
}

# Combine clinical and mutation matrices.
df <- merge(clinic, matrix %>% rownames_to_column("patient_id"), by = "patient_id") %>%
  column_to_rownames("patient_id")

df.combo <- df %>%
  dplyr::filter(treatment == "combo") %>%
  dplyr::select(-c(treatment:response)) %>%
  dplyr::filter(!is.na(cpr))

df.single <- df %>%
  dplyr::filter(treatment == "single") %>%
  dplyr::select(-c(treatment:response)) %>%
  dplyr::filter(!is.na(cpr))

# Keep genes with at least two mutation events in each cohort.
df.combo <- df.combo[, c(colnames(df.combo)[1], colnames(df.combo)[colSums(df.combo == "mut") >= 2])]
df.single <- df.single[, c(colnames(df.single)[1], colnames(df.single)[colSums(df.single == "mut") >= 2])]

# Run Fisher's exact test for each gene.
res.combo <- do.call(rbind, lapply(colnames(df.combo)[2:ncol(df.combo)], function(gene) {
  tab <- fisher.test(table(df.combo[, gene], df.combo[, "cpr"]))
  data.frame(gene = gene, odds = round(tab$estimate, 3), pval = round(tab$p.value, 4))
}))

res.single <- do.call(rbind, lapply(colnames(df.single)[2:ncol(df.single)], function(gene) {
  tab <- fisher.test(table(df.single[, gene], df.single[, "cpr"]))
  data.frame(gene = gene, odds = round(tab$estimate, 3), pval = round(tab$p.value, 4))
}))
