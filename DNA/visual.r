# visualize the genomic mutation landscape and mutation frequencies.

packages <- c("tidyverse", "openxlsx", "maftools")
invisible(lapply(packages, library, character.only = TRUE))

# Read clinical metadata.
sample <- read.csv("WES-metadata.csv", header = TRUE) %>%
  dplyr::filter(Tumor != "") %>%
  dplyr::filter(!grepl("p$", Patient))

meta <- read.xlsx("patient-info.xlsx", sheet = 1) %>%
  dplyr::filter(patient_id %in% sample$Patient) %>%
  mutate(
    treatment = ifelse(treatment == 2, "combo", "single"),
    orr = ifelse(orr == "ORR", "r", "nr"),
    group = paste(treatment, orr, sep = ".")
  ) %>%
  dplyr::select(patient_id, treatment, group)
colnames(meta) <- c("Tumor_Sample_Barcode", "treatment", "group")

# Read transcriptomic subtype annotations.
subtype <- read.csv("ccsubtype.csv", header = TRUE) %>%
  dplyr::select(-immune) %>%
  dplyr::filter(grepl("pre", sample_id)) %>%
  mutate(
    sample_id = str_sub(sample_id, 1, 5),
    subtype = case_when(
      subtype == 1 ~ "C2",
      subtype == 2 ~ "C4",
      subtype == 3 ~ "C1",
      subtype == 4 ~ "C3"
    )
  )
colnames(subtype) <- c("Tumor_Sample_Barcode", "subtype")

# Read mutation and genomic signature data.
load("vep.1000g.clean.rdata")
maf <- read.maf(maf.data, isTCGA = FALSE)

hrd <- read.table("scarhrd.txt", header = TRUE)
colnames(hrd) <- c("Tumor_Sample_Barcode", "HRD", "ALT", "LST", "HRD.sum")

msi <- read.table("msisensor2.txt", header = FALSE)
colnames(msi) <- c("Tumor_Sample_Barcode", "MSI")

# Merge clinical and genomic annotations.
meta <- meta %>%
  merge(subtype, by = "Tumor_Sample_Barcode", all.x = TRUE) %>%
  merge(hrd %>% dplyr::select(Tumor_Sample_Barcode, HRD.sum), by = "Tumor_Sample_Barcode", all.x = TRUE) %>%
  merge(msi, by = "Tumor_Sample_Barcode", all.x = TRUE)

# Define genes for the oncoplot.
target <- c(
  "TP53", "CDKN2A", "RB1", "CDK4",
  "KEAP1", "NOTCH1", "STK11",
  "ATM", "POLE", "BRCA2", "POLQ", "ATR",
  "KMT2D", "ARID1A", "KMT2C", "SMARCA4",
  "ERBB4", "EGFR", "ERBB2",
  "PIK3CA", "PTEN", "AKT1",
  "KRAS", "NF1", "BRAF"
)

# Build the mutation matrix for visualization.
length(unique(maf.data$Tumor_Sample_Barcode))
maf.data <- maf.data %>%
  dplyr::filter(!Tumor_Sample_Barcode %in% c("BT007", "BT016", "BT017", "BT047", "BT052", "BT115p", "BT122p"))
maf <- read.maf(maf.data, clinicalData = meta, isTCGA = FALSE)
str(maf@clinical.data)

mutation_cols <- c("#1F78B4", "#129490", "#70b77e", "#c5d86d", "#e3bd64", "#da9c6c", "#db6d6b", "#e07ba7")
names(mutation_cols) <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "In_Frame_Del",
  "Frame_Shift_Del", "In_Frame_Ins", "Frame_Shift_Ins", "Multi_Hit"
)

oncoplot(
  maf = maf,
  clinicalFeatures = c("group", "subtype", "HRD.sum", "MSI"),
  genes = target,
  colors = mutation_cols,
  annotationColor = list(
    group = c("combo.r" = "#f07167", "combo.nr" = "#fed9b7", "single.r" = "#00afb9", "single.nr" = "#0081a7"),
    subtype = c("C1" = "#bc4749", "C2" = "#a7c957", "C3" = "#6a994e", "C4" = "#386641"),
    HRD.sum = "Blues",
    MSI = "Reds"
  ),
  top = 100,
  fontSize = 0.5,
  keepGeneOrder = TRUE,
  sortByAnnotation = TRUE,
  removeNonMutated = FALSE
)

# Estimate mutation frequencies in single-agent and combination cohorts.
meta.single <- meta %>%
  dplyr::filter(treatment == "single")
barcode.single <- meta.single$Tumor_Sample_Barcode

res.single <- lapply(unique(maf@data$Hugo_Symbol), function(gene) {
  mutated_barcode <- unique(
    maf@data$Tumor_Sample_Barcode[
      maf@data$Hugo_Symbol == gene &
        maf@data$Tumor_Sample_Barcode %in% barcode.single
    ]
  )
  barcode <- maf@clinical.data %>%
    dplyr::filter(Tumor_Sample_Barcode %in% barcode.single) %>%
    pull(Tumor_Sample_Barcode)
  cluster <- maf@clinical.data %>%
    dplyr::filter(Tumor_Sample_Barcode %in% barcode.single) %>%
    pull(group)
  contingency <- table(factor(barcode %in% mutated_barcode, levels = c(FALSE, TRUE)), cluster)
  rownames(contingency) <- c("wild", "mutant")

  if (all(contingency >= 5)) {
    test <- chisq.test(contingency)
  } else {
    test <- fisher.test(contingency)
  }

  freq <- prop.table(contingency["mutant", ], margin = NULL) * 100
  cbind(data.frame(gene = gene, pval.single = test$p.value), t(freq))
}) %>%
  do.call(rbind, .)
res.single <- res.single %>%
  dplyr::filter(single.nr != "NaN" & single.r != "NaN")

meta.combo <- meta %>%
  dplyr::filter(treatment == "combo")
barcode.combo <- meta.combo$Tumor_Sample_Barcode

res.combo <- lapply(unique(maf@data$Hugo_Symbol), function(gene) {
  mutated_barcode <- unique(
    maf@data$Tumor_Sample_Barcode[
      maf@data$Hugo_Symbol == gene &
        maf@data$Tumor_Sample_Barcode %in% barcode.combo
    ]
  )
  barcode <- maf@clinical.data %>%
    dplyr::filter(Tumor_Sample_Barcode %in% barcode.combo) %>%
    pull(Tumor_Sample_Barcode)
  cluster <- maf@clinical.data %>%
    dplyr::filter(Tumor_Sample_Barcode %in% barcode.combo) %>%
    pull(group)
  contingency <- table(factor(barcode %in% mutated_barcode, levels = c(FALSE, TRUE)), cluster)
  rownames(contingency) <- c("wild", "mutant")

  if (all(contingency >= 5)) {
    test <- chisq.test(contingency)
  } else {
    test <- fisher.test(contingency)
  }

  freq <- prop.table(contingency["mutant", ], margin = NULL) * 100
  cbind(data.frame(gene = gene, pval.combo = test$p.value), t(freq))
}) %>%
  do.call(rbind, .)
res.combo <- res.combo %>%
  dplyr::filter(combo.nr != "NaN" & combo.r != "NaN")

# Combine mutation frequency results.
res <- merge(res.single, res.combo, by = "gene") %>%
  mutate(
    pval.single = as.numeric(pval.single),
    pval.combo = as.numeric(pval.combo),
    fdr.single = p.adjust(as.numeric(pval.single), method = "BH"),
    fdr.combo = p.adjust(as.numeric(pval.combo), method = "BH")
  ) %>%
  dplyr::filter(pval.single <= 0.1 | pval.combo <= 0.1)

# Inspect the mutation frequency for a representative gene.
gene.info <- getGeneSummary(maf) %>%
  dplyr::filter(Hugo_Symbol == "EPHA2")
gene.info$MutationFrequency <- gene.info$MutatedSamples / dplyr::n_distinct(maf@clinical.data$Tumor_Sample_Barcode)
