# calculate transcriptome-based immune signatures and compare response groups.

packages <- c("tidyverse", "openxlsx", "clusterProfiler", "GSVA", "ConsensusClusterPlus", "pheatmap", "ggplot2", "ggrepel")
invisible(lapply(packages, library, character.only = TRUE))

# Read TPM matrix and apply log2 transformation.
matrix <- read.csv("tpm.csv", header = TRUE, row.names = 1)
matrix <- log2(matrix + 1)

# Subset representative immune gene panels.
gene.immune <- matrix %>%
  dplyr::filter(rownames(.) %in% c(
    "PDCD1", "CD274", "CTLA4", "HAVCR2", "TIGIT", "LAG3", "BTLA", "TNFRSF14",
    "GZMB", "IFNG", "PRF1", "IL1B", "CXCL9", "CXCL10", "CXCL11", "CXCL13",
    "TRDC", "TRDV1", "TRDV2", "TRDV3", "KLRK1", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
    "B2M", "HLA-B", "HLA-C", "HLA-E", "TAP1", "CIITA"
  ))

gene.innate <- matrix %>%
  dplyr::filter(rownames(.) %in% c(
    "CGAS", "STING1", "TBK1", "IRF3", "IFNB1", "CXCL10",
    "TLR3", "TLR4", "TLR7", "TLR8", "TLR9", "TRIF",
    "RIGI", "MDA5", "MAVS",
    "NLRP3", "AIM2", "CASP1"
  ))

# Read immune deconvolution and subtype results.
immune <- read.csv("immune.csv", header = TRUE, row.names = 1)
subtype <- read.csv("ccsubtype.csv", header = TRUE)
cell.immune <- immune %>%
  dplyr::filter(grepl("tcellsi", rownames(.)))
colSums(cell.immune)

# Read TMB and clinical annotations.
load("tmb.rdata")
tmb <- tmb %>%
  mutate(sample_id = paste0(Tumor_Sample_Barcode, "T.pre")) %>%
  dplyr::select(sample_id, TMB)

clinic <- read.xlsx("sample-info.xlsx", sheet = 1) %>%
  dplyr::filter(sample_type %in% c("T"))

anno <- clinic %>%
  dplyr::select(sample_id, response, group1:group5) %>%
  merge(tmb, by = "sample_id") %>%
  merge(subtype %>% dplyr::select(sample_id, immune), by = "sample_id") %>%
  mutate(
    immune = case_when(
      immune == "Immune-Enriched, Non-Fibrotic" ~ "C1",
      immune == "Immune-Enriched, Fibrotic" ~ "C2",
      immune == "Fibrotic" ~ "C3",
      immune == "Depleted" ~ "C4"
    )
  ) %>%
  column_to_rownames("sample_id")

anno$group4 <- factor(
  anno$group4,
  levels = c(
    "combo.cpr.pre", "combo.cpr.post", "combo.ncpr.pre", "combo.ncpr.post",
    "single.cpr.pre", "single.cpr.post", "single.ncpr.pre", "single.ncpr.post"
  )
)

# Define immune-related gene signatures.
feature <- list(
  tls = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"),
  mhci.class = c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"),
  mhcii.class = c("CD74", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2"),
  tcell.activation = c("CXCL9", "CXCL10", "CXCL16", "IFNG", "IL15"),
  tcell.regulatory = c("FOXP3", "TNFRSF18"),
  tcell.exhaustion = c("CD274", "HAVCR2", "TNFRSF9", "CTLA4", "TOX"),
  tcell.survival = c("CD70", "CD27")
)
sig.feature <- gsva(as.matrix(matrix), feature, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

# Read public gene sets and calculate GSVA scores.
hallmark <- read.gmt("h.all.2024.Hs.gmt")
colnames(hallmark) <- c("term", "SYMBOL")
hallmark <- lapply(split(hallmark$SYMBOL, hallmark$term), function(x) x[x != ""])
sig.hallmark <- gsva(as.matrix(matrix), hallmark, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

ddr <- read.gmt("ddr.2024.Hs.gmt")
colnames(ddr) <- c("term", "SYMBOL")
ddr <- lapply(split(ddr$SYMBOL, ddr$term), function(x) x[x != ""])
sig.ddr <- gsva(as.matrix(matrix), ddr, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

apm <- read.gmt("apm.2024.Hs.gmt")
colnames(apm) <- c("term", "SYMBOL")
apm <- lapply(split(apm$SYMBOL, apm$term), function(x) x[x != ""])
sig.apm <- gsva(as.matrix(matrix), apm, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

innate <- read.gmt("immune.2024.Hs.gmt")
colnames(innate) <- c("term", "SYMBOL")
innate <- lapply(split(innate$SYMBOL, innate$term), function(x) x[x != ""])
sig.innate <- gsva(as.matrix(matrix), innate, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

nature <- read.gmt("nature.2024.Hs.gmt")
colnames(nature) <- c("term", "SYMBOL")
nature <- lapply(split(nature$SYMBOL, nature$term), function(x) x[x != ""])
sig.nature <- gsva(as.matrix(matrix), nature, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

# Scale signature scores for heatmap visualization.
table(colnames(sig.feature) == colnames(sig.apm))
df <- sig.innate %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame()

anno <- anno %>%
  arrange(group4, response, immune)
df <- df %>%
  dplyr::select(rownames(anno))
gap <- cumsum(table(anno$group4))[-length(table(anno$group4))]
bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))

pheatmap(
  df,
  annotation_col = anno,
  gaps_col = gap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  legend = FALSE,
  fontsize_row = 5,
  color = c(
    colorRampPalette(colors = c("#1e96fc", "#f8f9fa"))(length(bk) / 2),
    colorRampPalette(colors = c("#f8f9fa", "#f77f00"))(length(bk) / 2)
  ),
  annotation_colors = list(
    group4 = c(
      "combo.cpr.pre" = "#f2cc8f",
      "combo.cpr.post" = "#f07167",
      "combo.ncpr.pre" = "#ffd6a5",
      "combo.ncpr.post" = "#f28482",
      "single.cpr.pre" = "#80ed99",
      "single.cpr.post" = "#57cc99",
      "single.ncpr.pre" = "#48bfe3",
      "single.ncpr.post" = "#0081a7"
    ),
    response = c("cpr" = "#ed6a5a", "mpr" = "#ffd166", "non-mpr" = "#5ca4a9", "non-surgery" = "#adb5bd"),
    immune = c("C1" = "#bc4749", "C2" = "#a7c957", "C3" = "#6a994e", "C4" = "#386641")
  )
)

# Convert signature scores to long format for group comparison.
compare <- sig.nature %>%
  as.data.frame() %>%
  rownames_to_column("obj") %>%
  pivot_longer(-obj, names_to = "sample_id", values_to = "expr") %>%
  merge(clinic %>% dplyr::select(sample_id, group4), by = "sample_id")
colnames(compare) <- c("sample_id", "obj", "expr", "group")

compare <- compare %>%
  filter(!is.na(group))
compare$group <- factor(
  compare$group,
  levels = c(
    "combo.cpr.pre", "combo.cpr.post", "combo.ncpr.pre", "combo.ncpr.post",
    "single.cpr.pre", "single.cpr.post", "single.ncpr.pre", "single.ncpr.post"
  )
)

# Perform pairwise group comparisons for each signature.
res <- compare %>%
  group_by(obj) %>%
  summarise(
    combo.cpr_pval = wilcox.test(expr[group == levels(compare$group)[1]], expr[group == levels(compare$group)[2]], paired = FALSE, exact = FALSE)$p.value,
    combo.ncpr_pval = wilcox.test(expr[group == levels(compare$group)[3]], expr[group == levels(compare$group)[4]], paired = FALSE, exact = FALSE)$p.value,
    single.cpr_pval = wilcox.test(expr[group == levels(compare$group)[5]], expr[group == levels(compare$group)[6]], paired = FALSE, exact = FALSE)$p.value,
    single.ncpr_pval = wilcox.test(expr[group == levels(compare$group)[7]], expr[group == levels(compare$group)[8]], paired = FALSE, exact = FALSE)$p.value,
    combo.cpr_mean = log2(mean(expr[group == levels(compare$group)[2]]) / (mean(expr[group == levels(compare$group)[1]]) + 1e-6)),
    combo.ncpr_mean = log2(mean(expr[group == levels(compare$group)[4]]) / (mean(expr[group == levels(compare$group)[3]]) + 1e-6)),
    single.cpr_mean = log2(mean(expr[group == levels(compare$group)[6]]) / (mean(expr[group == levels(compare$group)[5]]) + 1e-6)),
    single.ncpr_mean = log2(mean(expr[group == levels(compare$group)[8]]) / (mean(expr[group == levels(compare$group)[7]]) + 1e-6)),
    .groups = "drop"
  ) %>%
  mutate(
    combo.cpr_padj = round(p.adjust(combo.cpr_pval, method = "BH"), 10),
    combo.ncpr_padj = round(p.adjust(combo.ncpr_pval, method = "BH"), 10),
    single.cpr_padj = round(p.adjust(single.cpr_pval, method = "BH"), 10),
    single.ncpr_padj = round(p.adjust(single.ncpr_pval, method = "BH"), 10)
  )

# Plot the bubble chart of differential signatures.
plot <- res %>%
  dplyr::select(obj, combo.cpr_mean:single.ncpr_mean, combo.cpr_padj:single.ncpr_padj)
colnames(plot) <- gsub("mean", "logfc", colnames(plot))

plot <- plot %>%
  pivot_longer(cols = -obj, names_to = c("comp", ".value"), names_sep = "_") %>%
  mutate(
    comp = factor(comp, levels = c("combo.cpr", "combo.ncpr", "single.cpr", "single.ncpr")),
    color = ifelse(logfc > 0, -log10(padj), -(-log10(padj))),
    size = abs(logfc)
  ) %>%
  dplyr::filter(!is.infinite(logfc) & !is.na(logfc) & logfc != "NaN")

summary(plot$size)
summary(plot$color)

ggplot(plot, aes(x = obj, y = comp)) +
  geom_point(aes(size = size, color = color), alpha = 0.9) +
  geom_text(data = subset(plot, padj < 0.1), aes(label = "*"), color = "black", size = 5, vjust = 0.7) +
  scale_size_continuous(name = "log2fc", range = c(2, 8), breaks = seq(0, 30, length.out = 5)) +
  scale_color_gradient2(name = "-log10(adj.P)", low = "#0096c7", mid = "white", high = "#d00000", midpoint = 0, limits = c(-4, 4)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 6),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.15),
    legend.position = "none"
  )
