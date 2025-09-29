# Preprocessing for Single-Cell Protein (CITE-seq) Data

This directory contains scripts to process raw CITE-seq data, converting it from single-cell resolution to pseudo-bulk format for DIVAS analysis.

## Prerequisites

- **Data Acquisition**: Run the data acquisition workflow to download raw CITE-seq files:
  ```bash
  cd ../../../data_acquisition
  bash download_arrayexpress_data.sh
  bash organize_arrayexpress_files.sh
  bash unzip_pro_data.sh
  ```
- **Data Extraction**: Extract `.txt` files from `data_acquisition/arrayexpress_data/pro_data_unzipped/`

## Processing Workflow

Execute the following scripts in order:

```bash
# 1. Create pseudo-bulk data matrix with CLR normalization
python create_pro_datablock.py
# Input:  ../../../data_acquisition/arrayexpress_data/pro_data_unzipped/*.txt
# Output: ../processed_omics_all/datablock_pro.tsv (CLR-normalized)

# 2. Filter for dual-timepoint cohort
python filter_sc_pro_for_120patients.py
# Input:  ../processed_omics_all/datablock_pro.tsv + sample_distribution/core_samples_*.csv
# Output: ../processed_omics_120/sc_pro_120patients.csv
```

## Key Parameters

Default processing settings:
- **Aggregation method**: Mean expression across cells (pseudo-bulk)
- **Normalization**: Centered Log Ratio (CLR) transformation (applied during pseudo-bulk creation)
- **Sample filtering**: Dual-timepoint patients only (T1 + T2)
- **Expected output**: ~25 proteins × 240 samples (120 patients × 2 timepoints)

## Integration

The final output (`sc_pro_120patients.csv`) is compatible with other omics datasets for DIVAS multi-omics integration analysis. The CLR normalization ensures proper handling of compositional protein expression data.
