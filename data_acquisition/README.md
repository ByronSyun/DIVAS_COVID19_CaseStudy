# Data Acquisition for COVID-19 Multi-Omics Study

This directory contains the scripts necessary to download and organize the raw single-cell sequencing data from the study publicly available on ArrayExpress under accession number [E-MTAB-11659](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-9357/sdrf). 

**Note**: The clinical metadata and bulk omics data (Proteomics, Metabolomics) are sourced from the associated Mendeley Data repository and should be downloaded manually. These scripts only handle the raw sequencing files from ArrayExpress (https://data.mendeley.com/datasets/tzydswhhb5/5).

## Step-by-Step Instructions

Execute the following three scripts in order from the `data_acquisition` directory:

```bash
# 1. Download raw single-cell data from ArrayExpress (E-MTAB-9357)
bash download_arrayexpress_data.sh

# 2. Organize files by data type
bash organize_arrayexpress_files.sh  

# 3. Unzip protein expression data
bash unzip_pro_data.sh
```

**Input:** ArrayExpress accession E-MTAB-9357 containing 1340 single-cell data files from Su et al. (2020) COVID-19 study.

**Output:** Organized directory structure:
```
arrayexpress_data/
├── gex_data/              # Single-cell gene expression
├── pro_data/              # CITE-seq protein expression  
├── pro_data_unzipped/     # Uncompressed protein data
├── cd8_tcr_data/          # CD8+ T-cell receptors
├── cd4_tcr_data/          # CD4+ T-cell receptors
└── bcr_data/              # B-cell receptors
```

The scripts will automatically create log files to track download and organization progress.
