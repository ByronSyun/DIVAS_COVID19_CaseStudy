# Multi-Omics Integration with DIVAS

This directory contains the 4-omics DIVAS integration analysis for COVID-19 data, combining single-cell RNA-seq, single-cell proteomics, bulk proteomics, and metabolomics.

## Overview

This analysis integrates four complementary omics datasets using DIVAS (Data Integration via Analysis of Subspaces) to identify joint and individual variation patterns across multiple biological scales.

## Input Data

The analysis uses preprocessed 120-patient datasets from `../preprocessing/processed_omics_120/`:
- **`sc_gex_120patients_aligned.csv`**: scRNA-seq data (8,634 genes × 240 samples)
- **`sc_pro_120patients.csv`**: sc-Proteomics data (25 proteins × 240 samples)  
- **`proteomics_120patients.csv`**: Bulk proteomics data (481 proteins × 240 samples)
- **`metabolomics_120patients.csv`**: Metabolomics data (763 metabolites × 240 samples)

## Quick Start

### Prerequisites
- **R** (>= 4.0) with DIVAS package
- **Required R packages**: DIVAS, CVXR, devtools

### Installation
```r
# Install DIVAS package
library(devtools)
install_github("ByronSyun/DIVAS_Develop/pkg", ref = "main", force = TRUE)
```

### Run Analysis
```bash
# Single command to run complete 4-omics DIVAS analysis
Rscript run_divas_analysis.R
```

## Analysis Details

### Input Data Specifications

| Omics Type | Features | Samples | Data Type | Source |
|------------|----------|---------|-----------|---------|
| **scRNA-seq** | 8,634 genes | 240 | Pseudo-bulk expression | Single-cell aggregation |
| **sc-Proteomics** | 25 proteins | 240 | CLR-normalized | CITE-seq antibodies |
| **Bulk Proteomics** | 481 proteins | 240 | Log-transformed | Plasma measurements |
| **Metabolomics** | 763 metabolites | 240 | Normalized | Plasma measurements |

### Sample Design
- **120 Patients**: Dual timepoint design (T1 + T2) = 240 samples total
- **Balanced Cohort**: Representative of COVID-19 severity spectrum
- **Timepoints**: T1 (baseline), T2 (follow-up)

### DIVAS Parameters
```r
divas_results <- DIVASmain(
  datablock = data_list,
  nsim = 100,           # Bootstrap simulations
  iprint = TRUE,        # Progress output
  colCent = TRUE,       # Column centering
  rowCent = FALSE,      # No row centering
  seed = 123,           # Reproducibility
  ReturnDetail = TRUE   # Detailed output
)
```

## Results Structure

### Main Results (`DIVAS_Results/divas_results_4omics_full.rds`)
```r
# Load results
results <- readRDS("DIVAS_Results/divas_results_4omics_full.rds")

# Key components:
results$Scores           # Sample scores matrix (240 × ~91 components)
results$Loadings         # Feature loadings per omics
results$LoadingsNames    # Joint component annotations
results$scoresList       # Individual omics scores
```

### Output Files
- **`divas_results_4omics_full.rds`**: Complete DIVAS output object
- **`divas_diagnostic_plots_4omics_full.pdf`**: Comprehensive diagnostic plots
- **`diagnostic_plots/`**: Individual diagnostic plots (PNG files)

## Expected Results
- **Total Components**: ~91 joint components identified
- **4-way Components**: Shared across all omics (highest biological relevance)
- **3-way Components**: Multi-omics disease signatures
- **2-way Components**: Pairwise omics relationships
- **Individual Components**: Omics-specific variation

## Citation

If you use this multi-omics integration analysis, please cite:

```
Prothero, J., ..., Marron J. S. (2024). 
Data integration via analysis of subspaces (DIVAS).
```
