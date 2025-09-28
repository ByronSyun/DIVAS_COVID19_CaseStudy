# CellTypist Cell Type Annotation

This directory contains the core scripts for performing cell type annotation using CellTypist PBMC model and preparing DIVAS input data.

## Overview

The workflow performs cell type annotation using CellTypist PBMC model. The annotated data is then directly used by the T1+T2 preparation workflow.

## Prerequisites

Before running the scripts in this directory, you must:

1. **Complete GEX preprocessing**: Run the single-cell GEX processing workflow to generate .h5ad files:
   ```bash
   cd ../../preprocessing/process_sc/sc_gex_processing
   python batch_process_gex.py
   ```

2. **Install required packages**:
   ```bash
   pip install scanpy celltypist pandas numpy tqdm
   ```

## Directory Structure

```
celltype_annotation/
├── run_celltypist_pbmc.py           # Main CellTypist annotation script
├── create_divas_input.py            # DIVAS input preparation script
├── Adult_COVID19_PBMC.pkl           # CellTypist PBMC model file
├── annotated_data_majority_voting/  # Output: Annotated .h5ad files (generated)
└── DIVAS_InputbyTime/              # Output: T1/T2 DIVAS input files (generated)
    ├── T1_B.csv
    ├── T1_CD14_Monocyte.csv
    ├── T1_CD4m_T.csv
    ├── T1_CD8m_T.csv
    ├── T1_NK.csv
    ├── T2_B.csv
    ├── T2_CD14_Monocyte.csv
    ├── T2_CD4m_T.csv
    ├── T2_CD8m_T.csv
    └── T2_NK.csv
```

## Step-by-Step Workflow

### Step 1: Run CellTypist PBMC Annotation

```bash
python run_celltypist_pbmc.py
```

**What it does:**
- Loads processed .h5ad files from `../../preprocessing/process_sc/sc_gex_processing/processed_gex_data`
- Applies CellTypist PBMC model for cell type annotation
- Uses majority voting to consolidate cell type predictions
- Saves annotated data to `annotated_data_majority_voting/`

**Expected output:**
- Multiple `*_annotated.h5ad` files in `annotated_data_majority_voting/`
- Annotation summary and statistics

### Step 2: Create DIVAS Input Files

```bash
python create_divas_input.py
```

**What it does:**
- Reads annotated .h5ad files from `annotated_data_majority_voting/`
- Identifies patients with both T1 and T2 timepoint data
- Extracts mean expression for 5 major cell types: CD14 Monocyte, B, CD8m T, CD4m T, NK
- Creates separate T1 and T2 data blocks for DIVAS analysis
- Saves data to `DIVAS_InputbyTime/`

**Expected output:**
- 10 CSV files (T1_*.csv and T2_*.csv for each cell type)
- Each file contains genes × samples matrix

## Cell Types Analyzed

The pipeline focuses on 5 major immune cell types:

1. **CD14 Monocyte** - Classical monocytes
2. **B** - B cells
3. **CD8m T** - CD8+ memory T cells  
4. **CD4m T** - CD4+ memory T cells
5. **NK** - Natural killer cells

## Data Flow

```
Raw .h5ad files (120 patients)
    ↓ [run_celltypist_pbmc.py]
Annotated .h5ad files
    ↓ [create_divas_input.py]
T1/T2 DIVAS input files (filtered patients with both timepoints)
    ↓ [Next: data_preparation/]
Combined 6-block DIVAS input
```

## Important Notes

### Large Data Files

- The `annotated_data_majority_voting/` directory will contain large .h5ad files (not included in git)
- These files are generated when you run the scripts and contain the full single-cell data
- Total size can be several GB depending on the dataset

### Sample Filtering

- Only patients with both T1 and T2 timepoint data are included in the final DIVAS input
- The exact number of patients may vary based on data availability and quality filters

### Memory Requirements

- CellTypist annotation can be memory-intensive for large datasets
- Monitor system memory usage during processing
- Consider processing samples in batches if memory is limited

## Troubleshooting

### Common Issues

1. **Missing input data**: Ensure GEX preprocessing is completed first
2. **Memory errors**: Reduce batch size or increase system memory
3. **Model download fails**: Check internet connection for CellTypist model download

### Expected Processing Time

- CellTypist annotation: 10-30 minutes per sample (depending on cell count)
- DIVAS input creation: 5-10 minutes total
- Total workflow: 1-3 hours for full dataset

## Next Steps

After completing this workflow, proceed to:

1. **Data Preparation**: `../data_preparation/` - Combine T1/T2 data for 6-block analysis
2. **DIVAS Analysis**: `../DIVAS_run/` - Run multi-omics integration analysis

## Output Files

The final output files in `DIVAS_InputbyTime/` are ready for downstream DIVAS analysis. Each CSV file contains:

- **Rows**: Gene features
- **Columns**: Patient samples 
- **Values**: Mean expression levels for the specific cell type and timepoint

Generated on: $(date)
