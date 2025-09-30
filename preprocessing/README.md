# Data Preprocessing Pipeline

This directory contains the complete data preprocessing pipeline for COVID-19 multi-omics data, preparing datasets for DIVAS analysis.

## Directory Structure & Data Flow

The preprocessing pipeline processes raw omics data through the following directories in order:

### 1. Single-Cell & Bulk Data Processing
- **`process_sc/`**: Processes single-cell data (scRNA-seq + CITE-seq) into pseudo-bulk format
- **`process_bulk/`**: Processes bulk omics data (proteomics + metabolomics)

### 2. Data Output Directories
- **`processed_omics_all/`**: Contains all processed samples (intermediate results)
- **`processed_omics_120/`**: Contains final 120-patient datasets (240 samples: T1 + T2 timepoints)
- **`sample_distribution/`**: Sample ID mapping and metadata files

## Processing Order

1. **Run single-cell processing**: `process_sc/sc_gex_processing/` and `process_sc/sc_pro_processing/`
2. **Run bulk processing**: `process_bulk/`
3. **Generated outputs**:
   - `processed_omics_all/`: All samples processed
   - `processed_omics_120/`: Final 120 dual-timepoint patients (240 samples) for DIVAS analysis

## Key Files

- **`sample_distribution/sample_ids.tsv`**: Sample ID conversion table
- **`sample_distribution/patients_meta.csv`**: Patient metadata and clinical information
- **`processed_omics_120/*.csv`**: Final aligned datasets ready for DIVAS integration

## Next Steps

After preprocessing completion:
1. **4-Omics Integration**: `../multi_omics_integration/run_divas_analysis.R`
2. **Cell Type Analysis**: `../scRNA_celltyist_analysis/`
