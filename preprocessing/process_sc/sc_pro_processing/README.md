# Preprocessing for Single-Cell Protein (CITE-seq) Data

This directory contains scripts to process the raw CITE-seq data, converting it from single-cell resolution to a pseudo-bulk format suitable for DIVAS analysis, and then filtering it for a specific patient cohort.

## Prerequisites

- You must first run the workflow in the `data_acquisition` directory to download and unzip the raw CITE-seq `.txt` files. These files should be located in `../../data_acquisition/arrayexpress_data/pro_data_unzipped/`.
- For Step 2, you must place the `core_samples.csv` file, which defines the patient cohort, into a subdirectory named `source_data` within this folder.

The required structure for Step 2 is:
```
sc_pro_processing/
├── source_data/
│   └── core_samples.csv
├── create_pro_datablock.py
└── filter_sc_pro_for_120patients.py
```

## Step-by-Step Workflow

Please execute the following scripts in order.

### Step 1: Create Pseudo-Bulk Data Matrix

This script aggregates the raw, single-cell protein expression data into a unified, sample-level matrix.

**Script:** `create_pro_datablock.py`

**What it Does:**
1.  **Reads Raw Data**: Scans the `pro_data_unzipped` directory for all individual sample `.txt` files.
2.  **Calculates Mean Expression**: For each sample, it computes the mean expression level for every protein across all cells, converting the data to a "pseudo-bulk" profile.
3.  **Standardizes Sample IDs**: Converts complex filenames into a clean, standardized format (e.g., `COVID_1_T1`, `Healthy_1053BW`).
4.  **Generates Output**:
    - `datablock_pro.tsv`: A comprehensive data matrix where rows are proteins and columns are all available samples.
    - `sample_ids.tsv`: A metadata file that maps the new standardized sample IDs back to their original filenames.

### Step 2: Filter for the 120-Patient Cohort

This script takes the full data matrix and filters it to include only the samples corresponding to the 120 patients with dual time points.

**Script:** `filter_sc_pro_for_120patients.py`

**What it Does:**
1.  **Reads Inputs**: Loads the full `datablock_pro.tsv`, the `sample_ids.tsv` mapping file, and the `core_samples.csv` list.
2.  **Matches Samples**: Identifies the exact sample columns in the datablock that correspond to the patients in the core list.
3.  **Filters Data**: Subsets the main data matrix, keeping only the columns for the matched samples.
4.  **Generates Output**:
    - `sc_pro_120patients_[timestamp].csv`: The final, filtered data matrix containing only the 120-patient cohort.
    - `sc_pro_120patients_metadata_[timestamp].csv`: The corresponding metadata for the filtered samples.
