# Data Preprocessing Pipeline

This directory contains the complete data preprocessing pipeline for COVID-19 multi-omics data, preparing datasets for DIVAS analysis.

## Overview

The preprocessing pipeline processes raw omics data and prepares standardized datasets for downstream DIVAS integration analysis. It handles four major data types:
- **Single-cell RNA-seq (scRNA-seq)**: Gene expression at single-cell resolution
- **Single-cell Proteomics (CITE-seq)**: Protein expression at single-cell resolution  
- **Bulk Proteomics**: Plasma protein measurements
- **Metabolomics**: Plasma metabolite measurements

## Directory Structure

```
preprocessing/
├── process_bulk/                    # Bulk omics data processing
│   ├── filter_metabolomics_for_120patients.py
│   ├── filter_proteomics_for_120patients.py
│   ├── improve_metabolomics_quality.R
│   ├── improve_proteomics_quality.R
│   ├── sample_ids.tsv              # Sample ID mapping (CRITICAL)
│   └── README.md
├── process_sc/                     # Single-cell data processing
│   ├── sc_gex_processing/          # scRNA-seq processing
│   │   ├── batch_process_gex.py
│   │   ├── filter_sc_gex_for_120patients.py
│   │   ├── gex_sample_metadata.tsv  # Sample metadata (CRITICAL)
│   │   └── README.md
│   ├── sc_pro_processing/          # sc-Proteomics processing
│   │   ├── filter_sc_pro_for_120patients.py
│   │   ├── normalize_sc_pro_clr.py
│   │   └── README.md
│   └── README.md
├── data_alignment/                 # Cross-omics sample alignment
│   ├── align_gex_final.py
│   ├── verify_final_alignment.py
│   └── README.md
├── sample_distribution/            # Sample metadata and selection
│   ├── core_samples_10X_metabolomics_proteomics_20250704_173841.csv
│   └── sample_ids.tsv
├── processed_omics_all/            # Intermediate results (all patients)
└── processed_omics_120/           # Final results (120 dual-timepoint patients)
    ├── metabolomics_120patients.csv
    ├── proteomics_120patients.csv
    ├── sc_gex_120patients_aligned.csv
    └── sc_pro_120patients.csv
```

## Processing Workflow

### Phase 1: Initial Data Processing (All Patients)
```bash
# Process bulk omics data
cd process_bulk
python process_metabolomics.py
python process_proteomics.py
R improve_metabolomics_quality.R
R improve_proteomics_quality.R

# Process single-cell data
cd ../process_sc/sc_gex_processing
python batch_process_gex.py

cd ../sc_pro_processing
python create_pro_datablock.py
```

### Phase 2: 120-Patient Filtering (Dual Timepoint)
```bash
# Filter for patients with both T1 and T2 data
cd process_bulk
python filter_metabolomics_for_120patients.py
python filter_proteomics_for_120patients.py

cd ../process_sc/sc_gex_processing  
python filter_sc_gex_for_120patients.py

cd ../sc_pro_processing
python filter_sc_pro_for_120patients.py
python normalize_sc_pro_clr.py
```

### Phase 3: Data Alignment
```bash
cd data_alignment
python align_gex_final.py
python verify_final_alignment.py
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

## Troubleshooting

### Common Issues

1. **Sample ID Mismatches**:
   ```bash
   # Check sample ID mapping
   head -5 sample_distribution/sample_ids.tsv
   ```

2. **Missing Data Files**:
   - Ensure raw data is placed in correct directories
   - Check file permissions and paths

3. **Memory Issues**:
   - Large datasets may require >16GB RAM
   - Consider processing in batches for very large cohorts

### Data Validation
```bash
cd data_alignment
python verify_final_alignment.py  # Validates sample consistency
```

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
