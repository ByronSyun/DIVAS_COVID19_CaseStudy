# Preprocessing for Bulk Proteomics and Metabolomics Data

This directory contains scripts to preprocess the bulk proteomics and metabolomics data for the COVID-19 study. The main goal is to extract the relevant data from the source files, clean them, and format them for the subsequent multi-omics integration analysis.

## Prerequisites

**1. Manual Data Download:**

Before running any scripts, you **must** manually download the supplementary data from the original publication's Mendeley Data repository:
- **Repository Link**: [Su, et al., Mendeley Data](https://data.mendeley.com/datasets/29ntw7g2wh/1)
- **File to Download**: `Table S1. Human subject details, plasma proteomic and metabolomic datasets and analysis...xlsx`

**2. Prepare Source Data:**

You need to extract the relevant sheets from the downloaded `.xlsx` file and save them as `.csv` files. Then, place them inside a `source_data` directory within this folder.

- Create a directory: `mkdir source_data`
- Prepare and place the following files inside `source_data/`:
  - `proteomics_data.csv`
  - `proteomics_metadata.csv`
  - `metabolomics_data.csv`
  - `metabolomics_metadata.csv`
  - `clinical_metadata.csv`
  - `healthy_donor_metadata.csv`

The final structure should look like this:
```
process_bulk/
├── source_data/
│   ├── proteomics_data.csv
│   ├── proteomics_metadata.csv
│   └── ... (other csv files)
├── run_data_processing.sh
├── process_proteomics.py
└── ... (other scripts)
```

## How to Run

Once the `source_data` directory is prepared, you can run the main processing script. This script will execute all the necessary Python and R scripts in the correct order.

```bash
bash run_data_processing.sh
```

## Detailed Workflow

The `run_data_processing.sh` script automates the following 5-step pipeline:

1.  **Step 1: Prepare Data (`update_covid_files.py`)**
    - Copies the necessary `.csv` files from the `source_data/` directory into the current working directory. This ensures the pipeline works on a local copy of the data.

2.  **Step 2: Process Metabolomics (`process_metabolomics.py`)**
    - Reads the `metabolomics_data.csv` and `metabolomics_metadata.csv` files.
    - Standardizes all sample IDs in both the data columns and the metadata `sample_id` column to a consistent format (e.g., from `INCOV001-BL` to `COVID_1_T1`).

3.  **Step 3: Process Proteomics (`process_proteomics.py`)**
    - Performs the same sample ID standardization for `proteomics_data.csv` and its associated metadata files.

4.  **Step 4: Quality Improvement for Metabolomics (`improve_metabolomics_quality.R`)**
    - This R script takes the processed metabolomics data and applies a rigorous quality control workflow:
    - **Missing Value Filtering**: Removes any metabolite (row) that has more than 20% missing values across all samples.
    - **KNN Imputation**: Uses the K-Nearest Neighbors algorithm (with k=10) to estimate and fill in any remaining missing values.
    - **PQN Normalization**: Applies Probabilistic Quotient Normalization to correct for sample-to-sample variations in sample dilution.
    - **PCA Visualization**: Generates a Principal Component Analysis (PCA) plot to provide a visual summary of the data quality and sample clustering after processing.

5.  **Step 5: Quality Improvement for Proteomics (`improve_proteomics_quality.R`)**
    - The same comprehensive quality improvement workflow (Filtering, KNN Imputation, PQN Normalization, and PCA) is applied to the proteomics data.

## Output

After the script finishes, two new subdirectories will be created:

-   **`improved_metabolomics/`**: Contains the final, high-quality metabolomics data (`improved_metabolomics_data.csv`), along with QC plots and reports.
-   **`improved_proteomics/`**: Contains the final, high-quality proteomics data (`improved_proteomics_data.csv`), along with QC plots and reports.

These directories hold the cleaned, analysis-ready data blocks for the subsequent multi-omics integration.
