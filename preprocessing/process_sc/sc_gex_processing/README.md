# Preprocessing for Single-Cell Gene Expression (scRNA-seq) Data

This directory contains scripts to process the raw scRNA-seq data, converting it from single-cell resolution to a pseudo-bulk format suitable for DIVAS analysis, and then filtering it for a specific patient cohort.

## Prerequisites

- You must first run the workflow in the `data_acquisition` directory to download and unzip the raw scRNA-seq `.tar.gz` files. These files should be extracted to `gex_data_unzipped/` directory.
- The `sample_ids.tsv` file must be available in the `../process_bulk/` directory for sample ID mapping.

The required structure for the complete workflow is:
```
sc_gex_processing/
├── gex_data_unzipped/          # Raw scRNA-seq data files (extracted from .tar.gz)
├── processed_gex_data/         # Output directory for processed samples
├── process_gex_data.py         # Main processing script
├── batch_process_gex.py        # Batch processing wrapper
├── create_divas_input.py       # DIVAS input preparation
├── check_data_quality.py      # Quality control checks
├── gex_sample_metadata.tsv    # Sample metadata mapping (generated)
├── gex_divas_datablock.tsv    # Final DIVAS input matrix (generated)
└── README.md                   # This documentation
```

## Step-by-Step Workflow

Please execute the following scripts in order.

### Step 1: Process Raw scRNA-seq Data

This script processes individual single-cell RNA-seq samples and converts them to pseudo-bulk format.

**Script:** `process_gex_data.py`

**What it Does:**
1. **Reads Raw Data**: Scans the `gex_data_unzipped` directory for all individual sample files.
2. **Quality Control**: Applies filtering based on gene expression thresholds and cell quality metrics.
3. **Pseudo-bulk Aggregation**: Calculates mean expression levels for each gene across all cells in each sample.
4. **Standardizes Sample IDs**: Converts complex filenames into standardized format (e.g., `COVID_1_T1`, `Healthy_1053BW`).
5. **Generates Output**:
   - Individual processed `.h5ad` files in `processed_gex_data/`
   - Quality control plots and statistics

### Step 2: Batch Processing (Optional)

For processing multiple samples efficiently, use the batch processing wrapper.

**Script:** `batch_process_gex.py`

**What it Does:**
- Automates the processing of all samples in the `gex_data_unzipped` directory
- Provides progress tracking and error handling
- Generates comprehensive processing logs

### Step 3: Quality Control Assessment

Check the quality of processed data before proceeding to integration.

**Script:** `check_data_quality.py`

**What it Does:**
1. **Data Integrity Checks**: Verifies that all expected samples were processed successfully
2. **Expression Statistics**: Calculates summary statistics for gene expression across samples
3. **Sample Coverage**: Checks coverage of genes and cells across the cohort
4. **Generates Reports**: Creates quality control reports and visualizations

### Step 4: Create DIVAS Input Files

This script prepares the final input files for DIVAS analysis.

**Script:** `create_divas_input.py`

**What it Does:**
1. **Aggregates Data**: Combines all processed samples into a unified expression matrix
2. **Gene Filtering**: Applies final filtering criteria for genes to be included in DIVAS
3. **Sample Matching**: Ensures sample IDs match with other omics data (proteomics, metabolomics)
4. **Generates Output**:
   - `gex_divas_datablock.tsv`: Main gene expression matrix for DIVAS
   - `gex_sample_metadata.tsv`: Sample metadata and mapping information

### Step 5: Filter for 120-Patient Cohort (Optional)

The `create_divas_input.py` script automatically handles sample ID mapping and filtering based on the sample mapping file. Additional filtering for specific patient cohorts can be done as a post-processing step if needed.

**Note:** The current scripts produce a comprehensive dataset with proper sample ID mapping. Manual filtering steps may be added later based on specific analysis requirements.

## File Descriptions

### Input Files
- **Raw scRNA-seq data**: Individual sample files in `gex_data_unzipped/`
- **sample_ids.tsv**: Sample ID mapping file (automatically located from `../process_bulk/`)

### Output Files
- **gex_divas_datablock.tsv**: Complete gene expression matrix (all samples)
- **sc_gex_120patients.tsv**: Filtered gene expression matrix (120-patient cohort)
- **gex_sample_metadata.tsv**: Sample metadata and ID mapping
- **processed_gex_data/**: Directory containing individual processed sample files
- **quality_reports/**: Directory containing data quality reports and visualizations

## Data Processing Parameters

The scripts use the following default parameters for quality control:
- **Minimum genes per cell**: 200
- **Minimum cells per gene**: 3
- **Maximum mitochondrial gene percentage**: 20%
- **Gene expression threshold**: Mean expression > 0.1 across samples

These parameters can be adjusted in the individual scripts as needed for your specific analysis requirements.

## Integration with Other Omics Data

The output from this workflow (`sc_gex_120patients.tsv`) is designed to be compatible with:
- Single-cell protein data (from `../sc_pro_processing/`)
- Bulk proteomics data
- Bulk metabolomics data

Ensure that sample IDs are consistent across all omics modalities before running DIVAS analysis.
