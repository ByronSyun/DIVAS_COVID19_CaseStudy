# DIVAS COVID-19 Case Study

📊 **[View the Interactive Case Study Website](https://byronsyun.github.io/DIVAS_COVID19_CaseStudy/)**

The interactive tutorial demonstrates DIVAS analysis using pre-computed results:
- Multi-omics data integration using DIVAS
- Component visualization and interpretation
- GO enrichment analysis
- Multi-omics correlation networks (Circos plots)

> **Note**: The tutorial focuses on DIVAS analysis and interpretation. For the complete data preprocessing pipeline (raw data download, quality control, cell type annotation, etc.), see the detailed instructions in the sections below.

---

## Data Availability

Due to file size limitations, large data files are hosted on Zenodo: https://doi.org/10.5281/zenodo.17430294

**Required download for tutorial:**
- `divas_results_combined_6block_renamed.rds` (419 MB) → place in `scRNA_celltyist_analysis/DIVAS_run/DIVAS_Results/`

**Optional:**
- `all_cells_metadata_complete.csv` (57.8 MB) → place in `scRNA_celltyist_analysis/celltype_annotation/`

All other data files (`Combined_*.csv`, `metadata.rds`) are already included in this repository.

---

## Overview

This repository contains a complete reproducible case study demonstrating the application of DIVAS (Data Integration via Analysis of Subspaces) to COVID-19 multi-omics data.

DIVAS is a novel statistical method for integrating multiple high-dimensional datasets. This case study applies DIVAS to COVID-19 patient data including:

- Single-cell RNA sequencing (scRNA-seq)
- Single-cell protein expression (CITE-seq) 
- Bulk proteomics
- Metabolomics

## Key Features

- **Complete Workflow**: From raw data processing to final DIVAS analysis
- **Reproducible**: All scripts with updated paths and dependencies
- **Well-documented**: Comprehensive README files in each directory
- **Cell Type Analysis**: Integration with CellTypist for automated annotation
- **Multi-timepoint**: Analysis of T1 and T2 timepoints (114 samples)

## Getting Started

### Prerequisites

- R (>= 4.0) with DIVAS package
- Python (>= 3.8) with scanpy, pandas, celltypist
- Required R packages: devtools, CVXR
- Required Python packages: scanpy, celltypist, pandas, numpy

### Installation

1. Clone this repository:
```bash
git clone https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy.git
cd DIVAS_COVID19_CaseStudy
```

2. Install DIVAS package:
```r
library(devtools)
install_github("ByronSyun/DIVAS_Develop/pkg", ref = "main")
```

3. Install Python dependencies:
```bash
pip install scanpy celltypist pandas numpy
```

### Usage

This case study follows a 4-phase workflow. Each phase has detailed instructions in its respective directory:

#### Phase 1: Data Acquisition
```bash
cd data_acquisition
# See data_acquisition/README.md for detailed instructions
bash download_arrayexpress_data.sh
```

#### Phase 2: Data Preprocessing  
```bash
cd preprocessing
# See preprocessing/README.md for complete pipeline
# Processes raw data → 120-patient standardized datasets
```

#### Phase 3: Multi-omics Integration (4-block DIVAS)
```bash
cd multi_omics_integration
# See multi_omics_integration/README.md for analysis details
Rscript run_divas_analysis.R
```

#### Phase 4: Cell Type Analysis (6-block DIVAS)
```bash
cd scRNA_celltyist_analysis
# See scRNA_celltyist_analysis/README.md for complete workflow
# Includes: CellTypist annotation → T1+T2 preparation → DIVAS analysis
```

> **Note**: Each directory contains a comprehensive `README.md` with detailed step-by-step instructions, troubleshooting guides, and expected outputs.

## Data Requirements for Full Preprocessing Pipeline

**Note**: The interactive tutorial (see link at top) uses pre-computed results from Zenodo and does **not** require running the preprocessing pipeline below.

If you want to reproduce the full preprocessing workflow from raw data:

**Required input data** (not included, must be downloaded separately):
- Raw scRNA-seq data (.h5ad format)
- Bulk proteomics data
- Metabolomics data  
- Sample metadata

**Generated data directories** (empty in Git, created during preprocessing):
- `preprocessing/processed_omics_120/`
- `preprocessing/processed_omics_all/`
- `scRNA_celltyist_analysis/celltype_annotation/annotated_data_majority_voting/`
- `multi_omics_integration/DIVAS_Results/`
- `scRNA_celltyist_analysis/DIVAS_run/divas_results/`

## Key Results

### 4-omics Integration
- Integrates scRNA-seq, sc-proteomics, bulk proteomics, and metabolomics
- 120 samples across two timepoints
- Identifies joint and individual variation patterns

### 6-block Cell Type Analysis  
- 114 samples (57 T1 + 57 T2)
- 4 cell types: CD4+ T, CD8+ T, CD14+ Monocytes, NK cells
- 2 bulk omics: proteomics, metabolomics
- 8,634 genes with comprehensive joint component analysis

## Citation

If you use this code or methodology, please cite:

```
Prothero, J., ..., Marron J. S. (2024). 
Data integration via analysis of subspaces (DIVAS).
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
