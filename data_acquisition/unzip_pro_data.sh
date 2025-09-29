#!/bin/bash

# Script to unzip all .txt.gz files in the pro_data directory

# Set paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
BASE_DIR="${SCRIPT_DIR}/arrayexpress_data"
SOURCE_DIR="${BASE_DIR}/pro_data"
TARGET_DIR="${BASE_DIR}/pro_data_unzipped"

# Make sure target directory exists
mkdir -p "${TARGET_DIR}"

echo "Starting to unzip pro_data files..."
echo "------------------------------------"
echo "Source directory: ${SOURCE_DIR}"
echo "Target directory: ${TARGET_DIR}"
echo "------------------------------------"

# Count total files
TOTAL_FILES=$(find "${SOURCE_DIR}" -name "*.txt.gz" | wc -l | tr -d ' ')
echo "Found ${TOTAL_FILES} files to unzip."
echo "------------------------------------"

# Initialize counter
COUNT=0

# Process each .txt.gz file
for gz_file in "${SOURCE_DIR}"/*.txt.gz; do
    # Extract the base filename without path and extension
    base_name=$(basename "${gz_file}" .gz)
    
    # Increment counter
    ((COUNT++))
    
    # Show progress
    echo "[${COUNT}/${TOTAL_FILES}] Unzipping: ${base_name}"
    
    # Unzip the file to the target directory
    gunzip -c "${gz_file}" > "${TARGET_DIR}/${base_name}"
done

echo "------------------------------------"
echo "Unzipping complete! Processed ${COUNT} files."
echo "Unzipped files are located in: ${TARGET_DIR}" 