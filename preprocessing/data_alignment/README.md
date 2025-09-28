# Data Alignment

This directory contains scripts for aligning and verifying the consistency of sample order across different omics data files.

## Overview

Before multi-omics integration analysis using DIVAS, it is crucial that all data matrices have identical sample ordering. This ensures proper alignment during analysis and prevents sample mismatches.

## Directory Structure

```
data_alignment/
├── align_gex_final.py          # Align GEX data to standard sample order
├── verify_final_alignment.py   # Verify sample order consistency
└── README.md                   # This file
```

## Scripts

### align_gex_final.py

**Purpose**: Aligns single-cell GEX data to match the sample order of other omics files.

**Input**: 
- Reference file: `../processed_omics_120/proteomics_120patients.csv` (for sample order)
- GEX data: `../processed_omics_all/gex_divas_datablock.tsv`

**Output**: 
- Aligned GEX file: `../processed_omics_120/sc_gex_120patients_aligned.csv`

**Usage**:
```bash
cd /Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/data_alignment
python align_gex_final.py
```

### verify_final_alignment.py

**Purpose**: Verifies that sample order is identical across all final omics data files.

**Input**: All files in `../processed_omics_120/`:
- `sc_gex_120patients_aligned.csv`
- `proteomics_120patients.csv`
- `metabolomics_120patients.csv`
- `sc_pro_120patients.csv`

**Output**: Console report showing verification results

**Usage**:
```bash
cd /Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/data_alignment
python verify_final_alignment.py
```

## Workflow

1. **After data filtering**: Run individual filtering scripts to generate 120-patient datasets
2. **Align GEX data**: Run `align_gex_final.py` to ensure GEX data matches sample order
3. **Verify alignment**: Run `verify_final_alignment.py` to confirm all files are properly aligned
4. **Proceed to DIVAS**: Only proceed with multi-omics analysis after successful verification

## Important Notes

- All omics data files must have identical sample ordering for DIVAS analysis
- The first column in all files should be named 'Feature' and contain feature IDs
- Sample IDs should follow the format: `COVID_X_TX` (e.g., `COVID_1_T1`, `COVID_1_T2`)
- If verification fails, re-run the alignment script before proceeding

## Dependencies

- pandas
- numpy (for data manipulation)

## Expected Output

After successful alignment and verification, you should see:
```
Success: All files have identical sample order
```

This confirms that all omics data files are ready for DIVAS multi-omics integration analysis.
