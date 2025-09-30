# Preprocessing for Single-Cell Gene Expression (scRNA-seq) Data

This directory contains scripts to process the raw scRNA-seq data, converting it from single-cell resolution to a pseudo-bulk format suitable for DIVAS analysis, and then filtering it for a specific patient cohort.

## Prerequisites

- **Data Acquisition**: Run the data acquisition workflow to download raw scRNA-seq files:
  ```bash
  cd ../../../data_acquisition
  bash download_arrayexpress_data.sh
  bash organize_arrayexpress_files.sh
  ```
- **Data Extraction**: Extract `.tar.gz` files from `data_acquisition/arrayexpress_data/gex_data/` to `gex_data_unzipped/` directory

## Processing Workflow

Execute the following scripts in order:

```bash
# 1. Filter raw data for dual-timepoint cohort (removes unwanted files)
python filter_sc_gex_for_120patients.py
# Input:  gex_data_unzipped/*.txt (all raw scRNA-seq files)
# Output: gex_data_unzipped/*.txt (filtered to 120-patient files only)

# 2. Process filtered scRNA-seq data → pseudo-bulk aggregation
# Option A: Process all files at once
python process_gex_data.py
# Input:  gex_data_unzipped/*.txt (120-patient filtered files)
# Output: processed_gex_data/*.h5ad (pseudo-bulk aggregated samples)

# Option B: Batch processing (recommended for large datasets)
python batch_process_gex.py  
# Input:  gex_data_unzipped/*.txt
# Output: processed_gex_data/*.h5ad (processes files in batches, calls process_gex_data.py)

# 3. Quality control assessment (Optional)
python check_data_quality.py
# Input:  processed_gex_data/*.h5ad
# Output: quality_reports/ (QC statistics and visualizations)

# 4. Create unified DIVAS input matrix
python create_divas_input.py
# Input:  processed_gex_data/*.h5ad + sample_distribution/sample_ids.tsv
# Output: gex_divas_datablock.tsv + gex_sample_metadata.tsv

# 5. Align GEX data with other omics datasets
python align_gex_final.py
# Input:  gex_divas_datablock.tsv + processed_omics_120/proteomics_120patients.csv
# Output: processed_omics_120/sc_gex_120patients_aligned.csv
```

## Key Parameters

Default quality control thresholds:
- **Minimum genes per cell**: 200 genes
- **Maximum genes per cell**: 2,500 genes
- **Minimum cells per gene**: 3 cells (≥0.1% of cells)
- **Maximum mitochondrial gene percentage**: 5%
- **Normalization target**: 10,000 counts per cell

## Integration

The complete workflow processes data through the following stages:

**Data Flow:**
- **`gex_data_unzipped/*.txt`**: Raw scRNA-seq files (filtered for 120-patient cohort)
- **`processed_gex_data/*.h5ad`**: Pseudo-bulk aggregated single-cell data (one .h5ad file per sample)
- **`gex_divas_datablock.tsv`**: 120-patient gene expression matrix for internal processing
- **`gex_sample_metadata.tsv`**: Sample metadata and ID mapping information
- **`../processed_omics_120/sc_gex_120patients_aligned.csv`**: Final aligned matrix compatible with other omics datasets (genes × 240 samples)

The alignment step (Step 5) ensures sample order consistency across all omics modalities, enabling seamless multi-omics integration analysis.
