# Data Preprocessing Pipeline

This directory contains the complete data preprocessing pipeline for COVID-19 multi-omics data, preparing datasets for DIVAS analysis.

## Overview

The preprocessing pipeline processes raw omics data and prepares standardized datasets for downstream DIVAS integration analysis. It focuses on **120 patients with dual timepoints (T1 baseline + T2 follow-up)** across four major data types:

- **Single-cell RNA-seq (scRNA-seq)**: Gene expression at single-cell resolution
- **Single-cell Proteomics (CITE-seq)**: Protein expression at single-cell resolution  
- **Bulk Proteomics**: Plasma protein measurements
- **Bulk Metabolomics**: Plasma metabolite measurements

The dual-timepoint design enables longitudinal analysis of COVID-19 disease progression and recovery patterns across multiple biological scales.

## Processing Workflow

Execute preprocessing in the following order:

### 1. Single-Cell Data Processing
- **scRNA-seq**: See `process_sc/sc_gex_processing/README.md` for detailed workflow
- **sc-Proteomics**: See `process_sc/sc_pro_processing/README.md` for detailed workflow

### 2. Bulk Omics Data Processing  
- **Bulk Proteomics & Metabolomics**: See `process_bulk/README.md` for detailed workflow

### 3. Final Data Verification
- **Sample Alignment**: See `processed_omics_120/README.md` for verification steps

## Directory Structure

```
preprocessing/
├── process_sc/                    # Single-cell data processing
│   ├── sc_gex_processing/        # scRNA-seq workflow
│   └── sc_pro_processing/        # sc-Proteomics workflow
├── process_bulk/                 # Bulk omics processing
├── processed_omics_120/          # Final 120-patient datasets
└── processed_omics_all/          # Intermediate results
```

## Key Features

- **120 dual-timepoint patients**: T1 baseline + T2 follow-up
- **4 omics modalities**: scRNA-seq, sc-Proteomics, Bulk Proteomics, Metabolomics  
- **Quality control**: Comprehensive data validation and filtering
- **Sample alignment**: Consistent ordering across all datasets
- **DIVAS-ready output**: Standardized format for multi-omics integration

## Next Steps

After successful preprocessing:
1. **4-Omics Integration**: `../multi_omics_integration/run_divas_analysis.R`
2. **Cell Type Analysis**: `../scRNA_celltyist_analysis/`

## Citation

This preprocessing pipeline was developed for:
```
Prothero, J., ..., Marron J. S. (2024). 
Data integration via analysis of subspaces (DIVAS).
```
