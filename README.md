# Code for LungMate-014

## Project Description

This repository contains the bioinformatics analyses supporting the translational component of the phase II clinical trial:

*Neoadjuvant toripalimab plus chemotherapy with or without tifcemalimab in potentially resectable stage III non-small cell lung cancer (LungMate-014): an open-label, single-center, randomized, phase II trial.*

The codebase focuses on downstream multi-omics analyses performed on trial-associated biospecimens, with emphasis on bulk transcriptomic profiling and DNA-based genomic characterization. The repository is organized as a script-based analytical workspace rather than a packaged end-to-end pipeline, and it documents the principal analyses underlying manuscript figures, integrative interpretation, and biomarker exploration.

## Repository Structure

### `RNA/`

The `RNA/` directory contains transcriptome-oriented analyses based primarily on bulk RNA-seq expression matrices and associated clinical/sample annotations. These scripts cover:

- gene identifier harmonization and duplicate-gene handling,
- count-to-TPM conversion and basic expression preprocessing,
- batch-related inspection and tumor purity estimation,
- sample organization and clinical feature integration,
- principal component analysis,
- transcriptomic subtype visualization,
- longitudinal subtype transition analysis,
- immune-feature comparisons across treatment and clinical strata.

### `DNA/`

The `DNA/` directory contains genomic analyses derived from mutation calls and trial-level sample annotations. These scripts cover:

- initial preprocessing and quality filtering of mutation data,
- driver-aware variant allele frequency filtering,
- integration of mutation data with clinical metadata,
- tumor mutational burden estimation,
- intratumor heterogeneity inference,
- gene-level association testing,
- mutation landscape visualization,
- integration of DNA-derived features with transcriptome-based immune signatures.

## Analysis Overview

### RNA-Based Analyses

The transcriptomic workflow begins with preprocessing of raw feature-level expression tables, including ENSEMBL-to-symbol conversion and removal of duplicated gene entries. Downstream scripts generate count and TPM matrices, align sample identifiers with clinical metadata, and support exploratory quality assessment.

Subsequent analyses focus on tumor transcriptomic characterization. These include principal component analysis, immune and molecular subtype annotation, GSVA/ssGSEA-based pathway scoring, subtype-specific heatmaps, and visualization of pre-treatment versus post-treatment shifts using bar plots and Sankey/alluvial diagrams. Additional scripts evaluate the balance of baseline clinical variables and compare immune-related expression features across treatment arms, PD-L1 strata, and genomic subgroups.

### DNA-Based Analyses

The genomic workflow starts from MAF-formatted mutation calls and associated sample metadata. Variants are filtered according to predefined mutation classes, PASS status, and variant allele frequency thresholds, with separate handling for known driver and non-driver genes. Cleaned mutation data are then used for downstream genomic summarization.

Derived genomic features include tumor mutational burden, mutation burden grouping, and intratumor heterogeneity metrics such as MATH-related summaries. These are integrated with clinical response annotations, HRD/MSI-derived features, and transcriptomic subtype information to support correlation analyses, group comparisons, and mutation landscape visualization. Additional scripts test gene-level mutation enrichment and examine DNA-associated immune and pathway features using transcriptome-derived signature scores.

## Methods Used in This Repository

The analytical methods implemented in the repository include standard statistical and visualization approaches commonly used in translational oncology studies, including:

- Fisher's exact test and chi-square test for categorical association analyses;
- Wilcoxon rank-sum test for non-parametric group comparisons;
- Spearman correlation for monotonic association analyses;
- principal component analysis for transcriptomic sample structure;
- GSVA/ssGSEA for pathway and immune-signature scoring;
- transcriptome-based immune deconvolution;
- `maftools`-based mutation summarization and oncoplot generation;
- heatmaps, bubble plots, boxplots, and Sankey/alluvial visualizations for integrative result presentation.

## Notes

- This repository is script-driven and intermediate-file-based; it is not implemented as a unified workflow engine.
- The code is intended to document and support the core bioinformatics analyses behind the translational findings of LungMate-014.
- Execution order is partly implicit and depends on the availability of intermediate objects produced by upstream scripts.
- Some scripts reflect manuscript-focused exploratory and figure-oriented analyses rather than a fully generalized software package.

## Context

This repository should be interpreted as the bioinformatics companion workspace for the translational analyses associated with LungMate-014. It is designed to support structured inspection of transcriptomic and genomic correlates of treatment exposure and response in potentially resectable stage III non-small cell lung cancer.
