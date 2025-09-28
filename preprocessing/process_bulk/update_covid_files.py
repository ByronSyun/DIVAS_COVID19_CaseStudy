#!/usr/bin/env python3
import os
import shutil

# This script prepares the workspace by copying the required source CSV files
# from the 'source_data' subdirectory into the current working directory.
# It is the first step in the preprocessing pipeline.

# Get the directory where the script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SOURCE_DIR = os.path.join(SCRIPT_DIR, "source_data")
TARGET_DIR = SCRIPT_DIR

print(f"Source data directory: {SOURCE_DIR}")
print(f"Target (output) directory: {TARGET_DIR}")

# Check if the source directory exists
if not os.path.isdir(SOURCE_DIR):
    print(f"Error: Source directory '{SOURCE_DIR}' not found.")
    print("Please create it and place the required CSV files inside, as described in the README.md.")
    exit(1)

# List of files to copy from the source directory
files_to_copy = [
    "metabolomics_data.csv",
    "metabolomics_metadata.csv",
    "proteomics_data.csv",
    "proteomics_metadata.csv",
    "clinical_metadata.csv",
    "healthy_donor_metadata.csv"
]

# Copy each file from the source to the target directory
copied_count = 0
for filename in files_to_copy:
    src_path = os.path.join(SOURCE_DIR, filename)
    dst_path = os.path.join(TARGET_DIR, filename)
    
    if os.path.exists(src_path):
        shutil.copy2(src_path, dst_path)
        copied_count += 1
    else:
        print(f"Warning: Source file not found, skipping: {src_path}")

if copied_count > 0:
    print(f"Successfully copied {copied_count} source files to the working directory.")
else:
    print("Error: No source files were found to copy. Please check the 'source_data' directory.")
    exit(1) 