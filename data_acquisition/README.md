# Data Acquisition for COVID-19 Multi-Omics Study

This directory contains the scripts necessary to download and organize the raw single-cell sequencing data from the study publicly available on ArrayExpress under accession number [E-MTAB-11659](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-9357/sdrf). 

**Note**: The clinical metadata and bulk omics data (Proteomics, Metabolomics) are sourced from the associated Mendeley Data repository and should be downloaded manually. These scripts only handle the raw sequencing files from ArrayExpress (https://data.mendeley.com/datasets/tzydswhhb5/5).

## Step-by-Step Instructions

Please execute the following scripts in order from the `DIVAS_COVID19_CaseStudy/data_acquisition` directory.

### Step 1: Download Raw Data

This script downloads all raw single-cell data files from the ArrayExpress FTP server.

**How to Run:**
```bash
bash download_arrayexpress_data.sh
```

**What it Does:**
- Reads a file list from `arrayexpress_data/processed-data_filelist.json` (Note: You may need to create this file or adjust the script if it's not present).
- Uses `wget` to download 1340 files corresponding to five data types:
  - Single-cell Gene Expression (GEX)
  - CITE-seq Protein Expression (PRO)
  - CD8+ T-cell Receptor (TCR)
  - CD4+ T-cell Receptor (TCR)
  - B-cell Receptor (BCR)
- Creates a log file `arrayexpress_data/download_arrayexpress.log`.

**Output:**
- A directory named `arrayexpress_data/` containing all 1340 raw data files (`.txt.gz`).

---

### Step 2: Unzip CITE-seq Protein Data

This script specifically unzips the CITE-seq (protein expression) data files, which are needed for downstream analysis.

**How to Run:**
```bash
bash unzip_pro_data.sh
```

**What it Does:**
- Finds all `.txt.gz` files in the `arrayexpress_data/pro_data/` directory (which will be created in the next step, so you may need to run Step 3 first or adjust paths).
- Unzips each file using `gunzip`.

**Output:**
- A new directory `arrayexpress_data/pro_data_unzipped/` containing the uncompressed `.txt` files for protein data.

---

### Step 3: Organize Downloaded Files

After downloading, this script organizes the 1340 raw files into separate subdirectories based on their data type.

**How to Run:**
```bash
bash organize_arrayexpress_files.sh
```

**What it Does:**
- Scans all `.txt.gz` files in the `arrayexpress_data/` directory.
- Moves each file into one of the following subdirectories based on its filename:
  - `gex_data/`
  - `pro_data/`
  - `cd8_tcr_data/`
  - `cd4_tcr_data/`
  - `bcr_data/`
- Creates a log file `arrayexpress_data/organize_files.log`.

**Output:**
- The `arrayexpress_data/` directory will now be cleanly organized with all files sorted into their respective data-type subfolders.
- **Final Structure:**
  ```
  arrayexpress_data/
  ├── gex_data/
  ├── pro_data/
  ├── cd8_tcr_data/
  ├── cd4_tcr_data/
  └── bcr_data/
  ```

After running these three scripts, all raw sequencing data will be downloaded and correctly structured for the next stages of the analysis.
