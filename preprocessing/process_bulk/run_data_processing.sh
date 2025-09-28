#!/bin/bash

# Main script to run the entire bulk data preprocessing workflow.
# This script executes Python and R scripts in sequence to clean,
# format, and perform quality control on the data.

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "=== Starting COVID-19 Bulk Data Preprocessing Workflow ==="
echo "Working Directory: ${SCRIPT_DIR}"
cd "${SCRIPT_DIR}"

# Step 1: Copy source files into the working directory
echo "--- Step 1: Preparing source data files... ---"
python3 update_covid_files.py
if [ $? -ne 0 ]; then
  echo "Error: Step 1 (update_covid_files.py) failed."
  exit 1
fi

# Step 2: Process metabolomics data
echo "--- Step 2: Processing metabolomics data... ---"
python3 process_metabolomics.py
if [ $? -ne 0 ]; then
  echo "Error: Step 2 (process_metabolomics.py) failed."
  exit 1
fi

# Step 3: Process proteomics data
echo "--- Step 3: Processing proteomics data... ---"
python3 process_proteomics.py
if [ $? -ne 0 ]; then
  echo "Error: Step 3 (process_proteomics.py) failed."
  exit 1
fi

# Step 4: Run quality control on metabolomics data
echo "--- Step 4: Running QC on metabolomics data... ---"
Rscript improve_metabolomics_quality.R
if [ $? -ne 0 ]; then
  echo "Error: Step 4 (improve_metabolomics_quality.R) failed."
  exit 1
fi

# Step 5: Run quality control on proteomics data
echo "--- Step 5: Running QC on proteomics data... ---"
Rscript improve_proteomics_quality.R
if [ $? -ne 0 ]; then
  echo "Error: Step 5 (improve_proteomics_quality.R) failed."
  exit 1
fi

echo "=== Bulk Data Preprocessing Workflow Completed Successfully ==="
echo "Final processed data is available in 'improved_metabolomics' and 'improved_proteomics' subdirectories." 