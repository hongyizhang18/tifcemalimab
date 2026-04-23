# perform initial filtering for genomic mutation data.

packages <- c("tidyverse", "openxlsx", "maftools")
invisible(lapply(packages, library, character.only = TRUE))

# Read driver gene annotations.
driver <- read.xlsx("gene-multiomics.xlsx", sheet = 5)
colnames(driver) <- driver[1, ]
driver <- driver[-1, ] %>%
  dplyr::filter(Cancer %in% c("LUAD", "LUSC"))

# Read clinical and sample metadata.
clinic <- read.xlsx("patient-info.xlsx", sheet = 1) %>%
  dplyr::select(patient_id, treatment, response)
sample <- read.csv("WES-metadata.csv", header = TRUE) %>%
  dplyr::filter(Tumor != "")

# Read MAF data and remap sample barcodes to patient IDs.
maf <- read.maf(maf = "vep.1000g.merge.maf")
maf.data <- maf@data
sample_map <- setNames(sample$Patient, sample$Tumor)

length(unique(maf.data$Tumor_Sample_Barcode))
maf.data$Tumor_Sample_Barcode <- maf.data$Tumor_Sample_Barcode %>%
  gsub("\\[|\\]", "", .) %>%
  ifelse(. %in% names(sample_map), sample_map[.], .)
length(unique(maf.data$Tumor_Sample_Barcode))

# Read RNA expression data for overlap checking.
rna <- read.csv("tpm.csv", header = TRUE, row.names = 1)
length(intersect(maf.data$Tumor_Sample_Barcode, unique(str_sub(colnames(rna), 1, 5))))

# Calculate variant allele frequency.
maf.data$t_vaf <- maf.data$t_alt_count / maf.data$t_depth

# Check supported mutation classes.
allowed_variant_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
  "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
  "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"
)
table(maf.data$Variant_Classification, useNA = "ifany")
if (!all(maf.data$Variant_Classification %in% allowed_variant_classes)) {
  stop("Unknown mutation types detected.")
}
cat("Primary:", nrow(maf.data), "rows remaining\n")

# Keep PASS variants only.
table(maf.data$FILTER, useNA = "ifany")
maf.data <- maf.data %>%
  dplyr::filter(FILTER == "PASS")
cat("FILTER:", nrow(maf.data), "rows remaining\n")

# Apply VAF thresholds for driver and non-driver genes.
summary(maf.data$t_vaf)
summary(maf.data$t_depth)
summary(maf.data$t_alt_count)
maf.data <- maf.data[
  is.na(maf.data$t_vaf) |
    (maf.data$Hugo_Symbol %in% driver$Gene & maf.data$t_vaf >= 0.015) |
    (!(maf.data$Hugo_Symbol %in% driver$Gene) & maf.data$t_vaf >= 0.05)
]
cat("VAF:", nrow(maf.data), "rows remaining\n")

# Save filtered mutation data.
save(maf.data, file = "vep.1000g.clean.rdata")
