# Multi-Omics Integration with DIVAS

This directory contains the 4-omics DIVAS integration analysis for COVID-19 data, combining single-cell RNA-seq, single-cell proteomics, bulk proteomics, and metabolomics.

## Overview

This analysis integrates four complementary omics datasets using DIVAS (Data Integration via Analysis of Subspaces) to identify:
- **Joint variation patterns** shared across multiple omics
- **Individual variation patterns** specific to each omics type
- **Multi-omics signatures** associated with COVID-19 progression

## Directory Structure

```
multi_omics_integration/
├── run_divas_analysis.R           # Main DIVAS analysis script
├── Final_Inputs/                  # Pre-processed 4-omics data (120 samples)
│   ├── sc_gex_120patients_aligned.csv      # scRNA-seq: 8,634 genes
│   ├── sc_pro_120patients.csv              # sc-Proteomics: 25 proteins  
│   ├── proteomics_120patients.csv          # Bulk Proteomics: 481 proteins
│   └── metabolomics_120patients.csv        # Metabolomics: 763 metabolites
├── DIVAS_Results/                 # Analysis outputs (generated)
│   ├── divas_results_4omics_full.rds       # Main DIVAS results
│   ├── divas_diagnostic_plots_4omics_full.pdf  # Diagnostic plots
│   └── diagnostic_plots/                   # Individual diagnostic plots
└── README.md                      # This file
```

## Quick Start

### Prerequisites
- **R** (>= 4.0) with DIVAS package
- **Required R packages**: DIVAS, CVXR, devtools
- **Preprocessed data**: 120-patient datasets in `Final_Inputs/`

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
| **scRNA-seq** | 8,634 genes | 120 | Pseudo-bulk expression | Single-cell aggregation |
| **sc-Proteomics** | 25 proteins | 120 | CLR-normalized | CITE-seq antibodies |
| **Bulk Proteomics** | 481 proteins | 120 | Log-transformed | Plasma measurements |
| **Metabolomics** | 763 metabolites | 120 | Normalized | Plasma measurements |

### Sample Design
- **120 Patients**: Dual timepoint design (T1 + T2)
- **Balanced Cohort**: Representative of COVID-19 severity spectrum
- **Timepoints**: T1 (baseline), T2 (follow-up)
- **Clinical Metadata**: WHO Ordinal Scale, demographics, outcomes

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

### Main Results (`divas_results_4omics_full.rds`)
```r
# Load results
results <- readRDS("DIVAS_Results/divas_results_4omics_full.rds")

# Key components:
results$Scores           # Sample scores matrix
results$Loadings         # Feature loadings per omics
results$LoadingsNames    # Joint component annotations
results$scoresList       # Individual omics scores
```

### Joint Components Analysis
- **Multi-omics components**: Shared variation across 2+ omics types
- **Individual components**: Omics-specific variation patterns
- **Component hierarchy**: From 4-way (all omics) to 1-way (individual)

### Diagnostic Outputs
- **Convergence plots**: Algorithm stability assessment
- **Variance explained**: Per-component contribution
- **Loading distributions**: Feature importance per omics

## Key Findings

### Expected Results
- **Total Components**: ~91 joint components identified
- **4-way Components**: Shared across all omics (highest biological relevance)
- **3-way Components**: Multi-omics disease signatures
- **2-way Components**: Pairwise omics relationships
- **Individual Components**: Omics-specific technical/biological variation

### Biological Interpretation
- **Immune Response Signatures**: Multi-omics patterns of immune activation
- **Metabolic Reprogramming**: Coordinated metabolic-proteomic changes
- **Cell Type Signatures**: scRNA-seq patterns reflected in bulk measurements
- **Disease Progression**: Temporal changes across timepoints

## Output Files

### Primary Results
- **`divas_results_4omics_full.rds`**: Complete DIVAS output object
- **`divas_diagnostic_plots_4omics_full.pdf`**: Comprehensive diagnostic plots

### Diagnostic Plots Directory
- **Individual PNG files**: Per-component diagnostic plots
- **Convergence plots**: Algorithm performance assessment
- **Loading plots**: Feature contribution visualization

## Usage Examples

### Loading and Exploring Results
```r
# Load results
library(DIVAS)
results <- readRDS("DIVAS_Results/divas_results_4omics_full.rds")

# Examine structure
names(results)
dim(results$Scores)  # Sample scores: 120 × 91

# Check joint components
head(results$LoadingsNames$sc_gex)  # scRNA-seq component names

# Extract 4-way components (shared across all omics)
four_way_components <- grep("sc_gex\\+sc_pro\\+proteomics\\+metabolomics", 
                           results$LoadingsNames$sc_gex, value = TRUE)
```

### Visualization
```r
# Plot sample scores for first component
plot(results$Scores[,1], results$Scores[,2], 
     main = "DIVAS Sample Scores", 
     xlab = "Component 1", ylab = "Component 2")

# Examine loadings for a specific omics
heatmap(results$Loadings$sc_gex[,1:10], 
        main = "scRNA-seq Loadings")
```

## Performance Metrics

### Computational Requirements
- **Runtime**: ~10-15 minutes (depends on system)
- **Memory**: ~8-16GB RAM recommended
- **Storage**: ~500MB for complete results

### Quality Indicators
- **Convergence**: All components should converge (check diagnostic plots)
- **Variance Explained**: Components should capture meaningful variation
- **Biological Relevance**: Multi-omics components should show coherent patterns

## Troubleshooting

### Common Issues

1. **Memory Errors**:
   ```r
   # Reduce bootstrap simulations if needed
   nsim = 50  # instead of 100
   ```

2. **Package Installation Issues**:
   ```r
   # Install dependencies first
   install.packages(c("CVXR", "devtools"))
   # Then install DIVAS
   devtools::install_github("ByronSyun/DIVAS_Develop/pkg", ref = "main", force = TRUE)
   ```

3. **Data Loading Issues**:
   - Ensure all 4 CSV files are present in `Final_Inputs/`
   - Check file permissions and paths
   - Verify data format (samples as columns, features as rows)

### Validation Checks
```r
# Verify input data consistency
data_files <- list.files("Final_Inputs", pattern = "*.csv", full.names = TRUE)
for(file in data_files) {
  data <- read.csv(file, row.names = 1)
  cat(basename(file), ": ", nrow(data), "×", ncol(data), "\n")
}
```

## Integration with Other Analyses

### Upstream Dependencies
- **Preprocessing Pipeline**: `../preprocessing/` must be completed first
- **Sample Alignment**: Requires consistent sample ordering across omics

### Downstream Applications
- **Cell Type Analysis**: Results inform `../scRNA_celltyist_analysis/`
- **Pathway Analysis**: Joint components can be analyzed for biological pathways
- **Clinical Associations**: Sample scores can be correlated with clinical outcomes

## Citation

If you use this multi-omics integration analysis, please cite:

```
Prothero, J., ..., Marron J. S. (2024). 
Data integration via analysis of subspaces (DIVAS).
```

## Contact

For questions about this analysis:
- **Issues**: https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy/issues
- **DIVAS Package**: https://github.com/ByronSyun/DIVAS_Develop
