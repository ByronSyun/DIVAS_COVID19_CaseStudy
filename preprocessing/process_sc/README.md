# Single-Cell Data Preprocessing

This directory contains workflows for preprocessing single-cell omics data for DIVAS analysis. The preprocessing converts single-cell resolution data to pseudo-bulk format suitable for multi-omics integration.

## Directory Structure

```
process_sc/
├── sc_gex_processing/          # Single-cell RNA-seq data preprocessing
│   ├── process_gex_data.py     # Main GEX processing script
│   ├── batch_process_gex.py    # Batch processing wrapper
│   ├── create_divas_input.py   # DIVAS input preparation
│   ├── check_data_quality.py  # Quality control checks
│   ├── gex_data_unzipped/      # Raw scRNA-seq data (user provides)
│   ├── processed_gex_data/     # Processed output directory
│   ├── gex_sample_metadata.tsv # Sample metadata
│   ├── sc_gex_120patients.tsv  # Final filtered dataset
│   └── README.md               # Detailed GEX workflow documentation
│
└── sc_pro_processing/          # Single-cell protein (CITE-seq) preprocessing
    ├── create_pro_datablock.py # Main protein processing script
    ├── filter_sc_pro_for_120patients.py # Patient cohort filtering
    └── README.md               # Detailed protein workflow documentation
```

## Workflow Overview

### Phase 1: Single-Cell Protein Processing
1. **Raw Data**: CITE-seq `.txt` files from ArrayExpress
2. **Processing**: Aggregate protein expression by sample (pseudo-bulk)
3. **Filtering**: Filter for 120-patient cohort
4. **Output**: `sc_pro_120patients_[timestamp].csv`

### Phase 2: Single-Cell RNA Processing  
1. **Raw Data**: scRNA-seq `.tar.gz` files from ArrayExpress
2. **Processing**: Quality control, normalization, pseudo-bulk aggregation
3. **Integration**: Create unified gene expression matrix
4. **Filtering**: Filter for 120-patient cohort
5. **Output**: `sc_gex_120patients.tsv`

## Prerequisites

### Data Requirements
Both workflows require raw data from ArrayExpress. Please first complete the data acquisition workflow:

```bash
# From the root directory
cd DIVAS_COVID19_CaseStudy/data_acquisition/
./download_arrayexpress_data.sh
./organize_arrayexpress_files.sh
```

### Patient Cohort Definition
Both workflows require the `core_samples.csv` file that defines the 120-patient cohort. This file should be placed in a `source_data/` subdirectory within each processing folder.

## Execution Order

1. **Complete data acquisition** (see `data_acquisition/README.md`)
2. **Process single-cell proteins** (see `sc_pro_processing/README.md`)
3. **Process single-cell RNA** (see `sc_gex_processing/README.md`)

## Configuration

All Python scripts have been configured with the correct paths for the new directory structure. The scripts are ready to run without additional path modifications.

## Integration with DIVAS Analysis

The outputs from both workflows are designed to be compatible with the bulk omics data and can be used together in the multi-omics integration phase:

- **Bulk omics**: Proteomics, metabolomics (from `../process_bulk/`)
- **Single-cell omics**: scRNA-seq, CITE-seq proteins (from this directory)
- **Integration**: Combined analysis using DIVAS methodology

## Quality Control

Both workflows include quality control steps:
- Data completeness checks
- Expression statistics and distributions  
- Sample coverage analysis
- Integration compatibility verification

Refer to the individual README files in each subdirectory for detailed workflow documentation.
