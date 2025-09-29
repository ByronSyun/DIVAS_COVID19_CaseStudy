# Data Preprocessing Pipeline

This directory contains the complete data preprocessing pipeline for COVID-19 multi-omics data, preparing datasets for DIVAS analysis.

## Overview

The preprocessing pipeline processes raw omics data and prepares standardized datasets for downstream DIVAS integration analysis. It handles four major data types:
- **Single-cell RNA-seq (scRNA-seq)**: Gene expression at single-cell resolution
- **Single-cell Proteomics (CITE-seq)**: Protein expression at single-cell resolution  
- **Bulk Proteomics**: Plasma protein measurements
- **Bulk Metabolomics**: Plasma metabolite measurements

## Processing Workflow

### Single-Cell Data Processing
```bash
# 1. Process scRNA-seq data
cd process_sc/sc_gex_processing
python batch_process_gex.py                    # Process raw 10X data
python filter_sc_gex_for_120patients.py        # Filter for dual-timepoint patients

# 2. Process single-cell proteomics (CITE-seq)
cd ../sc_pro_processing
python create_pro_datablock.py                 # Create protein datablock
python filter_sc_pro_for_120patients.py        # Filter for dual-timepoint patients
python normalize_sc_pro_clr.py                 # Apply CLR normalization
```

### Bulk Omics Data Processing
```bash
# 1. Process bulk proteomics and metabolomics
cd process_bulk
python process_metabolomics.py                 # Initial metabolomics processing
python process_proteomics.py                   # Initial proteomics processing
Rscript improve_metabolomics_quality.R         # Quality improvement
Rscript improve_proteomics_quality.R           # Quality improvement

# 2. Filter for dual-timepoint patients
python filter_metabolomics_for_120patients.py  # Extract 120 patients
python filter_proteomics_for_120patients.py    # Extract 120 patients
```

### Data Alignment and Finalization
```bash
# Ensure sample order consistency across all omics
cd data_alignment
python align_gex_final.py                      # Align scRNA-seq sample order
python verify_final_alignment.py               # Verify consistency
```

## Key Features

### 🔄 **Dual Processing Strategy**
- **All Patients**: Initial processing of complete dataset
- **120 Patients**: Focused analysis on patients with dual timepoints (T1 + T2)

### 📊 **Quality Control**
- **Data Validation**: Comprehensive checks for missing values, outliers
- **Sample Alignment**: Ensures consistent sample ordering across omics
- **ID Mapping**: Robust sample ID standardization across datasets

### 🎯 **DIVAS-Ready Output**
- **Standardized Format**: All outputs in CSV format with consistent structure
- **Sample Consistency**: Identical sample ordering across all omics types
- **Feature Selection**: Quality-filtered features ready for integration

## Data Requirements

### Input Data (Required)
- **Raw scRNA-seq data**: 10X Genomics format or H5AD files
- **Raw proteomics data**: Protein abundance matrices
- **Raw metabolomics data**: Metabolite abundance matrices
- **Sample metadata**: Patient demographics, timepoints, clinical data

### Output Data (Generated)
- **processed_omics_120/**: Final 120-patient datasets ready for DIVAS
- **processed_omics_all/**: Intermediate results for all patients

## Sample Selection Criteria

### 120-Patient Cohort
- **Dual Timepoint**: Patients with both T1 (baseline) and T2 (follow-up) samples
- **Complete Data**: Availability across all 4 omics modalities
- **Quality Passed**: Samples passing all QC filters

### Sample ID Standardization
- **Consistent Format**: `COVID_XXX_T1` / `COVID_XXX_T2`
- **Cross-Reference**: `sample_ids.tsv` provides mapping between formats
- **Validation**: Automated checks ensure ID consistency

## Key Files

### Critical Configuration Files
- **`sample_ids.tsv`**: Master sample ID mapping (DO NOT MODIFY)
- **`gex_sample_metadata.tsv`**: scRNA-seq sample metadata
- **`core_samples_*.csv`**: Patient selection criteria

### Final Output Files
- **`*_120patients.csv`**: Ready for DIVAS analysis
- **Data aligned**: Samples in identical order across all files


## Quality Metrics

### Expected Output Sizes
- **Metabolomics**: ~763 metabolites × 120 samples
- **Proteomics**: ~481 proteins × 120 samples  
- **scRNA-seq**: ~8,634 genes × 120 samples
- **sc-Proteomics**: ~25 proteins × 120 samples

### Success Criteria
- ✅ All 4 omics files have identical sample ordering
- ✅ No missing values in critical samples
- ✅ Sample IDs match across all datasets
- ✅ Feature counts within expected ranges

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
