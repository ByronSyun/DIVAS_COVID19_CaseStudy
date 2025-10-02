# Single-cell RNA-seq Cell Type Analysis with DIVAS

This directory contains the complete workflow for single-cell RNA-seq cell type analysis using CellTypist annotation and DIVAS integration for 114 balanced samples.

## Overview

This analysis is based on the multi-omics COVID-19 study by Su et al. (2020), which identified distinct immune cell subpopulations associated with disease severity. Rather than analyzing the entire single-cell GEX data as a single block, we extract specific cell types to understand their individual contributions to the immune response.

**Analysis workflow**:
1. **Annotates cell types** using CellTypist PBMC model
2. **Combines T1+T2 data** for 4 major cell types (114 samples, ~8.5k common genes)
3. **Prepares 6-block DIVAS input** (4 cell types: ~8.5k genes each + 2 bulk omics: 481 proteins + 763 metabolites)
4. **Ready for DIVAS analysis** to identify coordinated patterns across immune cell types and molecular data

**Reference**: Su, Y., Chen, D., Yuan, D., et al. (2020). Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19. *Cell*, 183(6), 1479-1495.

## Complete Workflow

### Step 1: Cell Type Annotation

**Input Data**: Processed single-cell GEX data from `../preprocessing/process_sc/sc_gex_processing/processed_gex_data/`
- **Format**: Individual `.h5ad` files for each sample (120 patients × 2 timepoints)
- **Content**: Quality-filtered single-cell gene expression data

**CellTypist PBMC Model**:
- **Model**: Adult COVID-19 PBMC reference model (`Adult_COVID19_PBMC.pkl`)
- **Framework**: Logistic regression-based automated cell type annotation
- **Reference**: Domínguez Conde et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*, 376(6594), eabl5197.
- **Documentation**: [CellTypist Models](https://www.celltypist.org/models)
- **Method**: Machine learning-based cell type classification with majority voting for robust annotations

**Run annotation**:
```bash
cd celltype_annotation
python run_celltypist_pbmc.py
```

**Output**: Cell type-annotated `.h5ad` files in `celltype_annotation/annotated_data_majority_voting/`
- **Format**: AnnData objects with cell type labels in `.obs['majority_voting']`
- **Size**: Several GB (not included in git repository)

### Step 2: T1+T2 Data Combination

**Sample Selection**: The analysis uses 114 carefully selected samples from `T1T2_preparation/balanced_sample_selection_114.csv`:

- **Selection Strategy**: Balanced sampling across COVID-19 severity and clinical metadata
- **WHO Ordinal Scale Balance**: Proportional representation across WOS scores (0-7)
- **Severity Balance**: Equal representation of healthy, mild, moderate, and severe cases
- **Temporal Design**: 57 T1 (baseline) + 57 T2 (follow-up) samples
- **Clinical Metadata**: Includes sex, age, ethnicity for comprehensive analysis

*Note: The sample selection algorithm was performed during the original analysis phase to ensure optimal balance for DIVAS correlation sensitivity.*

**Target Cell Types** (5 major immune populations):
1. **CD14 Monocyte** - Classical monocytes
2. **B** - B cells  
3. **CD8 T** - CD8+ T cells (includes CD8m, CD8n subtypes)
4. **CD4 T** - CD4+ T cells (includes CD4m, CD4n subtypes)
5. **NK** - Natural killer cells

**Data Combination Process**:
```bash
cd T1T2_preparation
Rscript combine_t1_t2_data.R
```

**Key Processing Steps**:
1. **Sample Selection**: Uses `balanced_sample_selection_114.csv` to identify target samples
2. **Cell Type Aggregation**: Combines specific subtypes into broader categories (e.g., CD4m T + CD4n T → CD4 T)
3. **Gene Filtering**: Applies common gene intersection and 10% expression threshold per timepoint
4. **T1/T2 Integration**: Combines T1 and T2 data while maintaining temporal information
5. **Bulk Data Integration**: Incorporates proteomics and metabolomics from `../preprocessing/processed_omics_120/`

**Output**: 6-block DIVAS input files in `DIVAS_run/divas_input_T1T2combined/`:
- **`Combined_CD14_Monocyte.csv`**: CD14+ monocyte expression (8634 genes × 114 samples)
- **`Combined_CD4_T_combined.csv`**: CD4+ T cell expression (8634 genes × 114 samples)
- **`Combined_CD8_T_combined.csv`**: CD8+ T cell expression (8634 genes × 114 samples)
- **`Combined_NK.csv`**: NK cell expression (8634 genes × 114 samples)
- **`Combined_proteomics.csv`**: Bulk proteomics (481 proteins × 114 samples)
- **`Combined_metabolomics.csv`**: Bulk metabolomics (763 metabolites × 114 samples)

### Data Summary

| Data Block | Features | Total Samples | T1 Samples | T2 Samples |
|------------|----------|---------------|------------|------------|
| **CD4_T_combined** | 8,634 | 114 | 57 | 57 |
| **CD8_T_combined** | 8,634 | 114 | 57 | 57 |
| **CD14_Monocyte** | 8,634 | 114 | 57 | 57 |
| **NK** | 8,634 | 114 | 57 | 57 |
| **proteomics** | 481 | 114 | 57 | 57 |
| **metabolomics** | 763 | 114 | 57 | 57 |

*Perfect temporal balance: Each data block contains exactly 57 T1 (baseline) and 57 T2 (follow-up) samples for robust longitudinal analysis.*

## Key Features

### 🔬 **Cell Type Analysis**
- **Automated Annotation**: CellTypist PBMC model with majority voting
- **Quality Control**: Sample coverage and cell type distribution checks
- **Broad Categories**: Aggregates specific subtypes for robust analysis

### 📊 **6-block DIVAS Integration**
- **scRNA-seq Blocks**: 4 cell type-specific expression profiles (8634 genes each)
- **Bulk Omics Blocks**: Proteomics (481 features) + Metabolomics (763 features)
- **Sample Design**: 114 samples from balanced clinical selection
- **Gene Coverage**: 8634 genes after common gene filtering and quality control

### ⚡ **Streamlined Workflow**
- **Direct Processing**: Reads directly from `.h5ad` files
- **Pre-computed Selection**: Uses optimized 114-sample configuration
- **Ready for DIVAS**: Output files formatted for immediate DIVAS analysis

## Data Requirements

### Prerequisites
- **R** (>= 4.0) with required packages (dplyr)
- **Python** (>= 3.8) with scanpy, celltypist, pandas
- **Preprocessed Data**: 
  - `.h5ad` files from GEX preprocessing
  - Bulk omics data from `../preprocessing/processed_omics_120/`

### Large Files (Not in Git)
- **`annotated_data_majority_voting/`**: Several GB of annotated single-cell data
- **Input `.h5ad` files**: Must be generated by preprocessing workflow

### Included Files
- **`balanced_sample_selection_114.csv`**: Pre-computed sample selection (28KB)
- **`divas_input_T1T2combined/`**: Final 6-block data (44MB, generated by workflow)

## Important Notes

### Sample Selection
The 114 samples were selected using a sophisticated balancing algorithm that:
- Ensures representation across COVID-19 severity spectrum
- Maintains WHO Ordinal Scale distribution balance
- Accounts for temporal pairing (T1/T2 availability)
- Optimizes for DIVAS correlation sensitivity analysis

### Gene Filtering
The final 8634 genes represent:
- **Common genes**: Present across all 114 samples
- **Expression threshold**: ≥10% expression within each timepoint
- **T1/T2 intersection**: Genes meeting criteria in both timepoints
- **Quality control**: Removes low-quality and sparse features

## Metadata Files
- **`final_sample_info.csv`**: Sample metadata with clinical information
- **`data_summary.csv`**: Data block specifications and dimensions

### Memory Requirements
- **CellTypist annotation**: 8-16GB RAM recommended
- **Data combination**: 16-32GB RAM for large `.h5ad` processing
- **Storage**: ~50GB for intermediate files (not in git)

## Next Steps

After completing this workflow:
1. **DIVAS Analysis**: Use `DIVAS_run/run_combined_6block_divas.R` (to be organized as R Markdown)
2. **Results Visualization**: Generate diagnostic plots and component analysis
3. **Biological Interpretation**: Analyze joint components and multi-omics patterns

## Citation

If you use this cell type analysis workflow, please cite:

```
Prothero, J., ..., Marron J. S. (2024). 
Data integration via analysis of subspaces (DIVAS).
```

---

*This workflow represents the core single-cell analysis pipeline for the COVID-19 multi-omics DIVAS study, focusing on 114 balanced samples for optimal statistical power and biological interpretation.*