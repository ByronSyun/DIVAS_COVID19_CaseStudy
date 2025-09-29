#!/bin/bash

# Script to download single-cell data from ArrayExpress for E-MTAB-9357
# This script reads a JSON file listing all files to download.

# Ensure the script is run from the workspace root or adjust paths accordingly
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # Should be DIVAS_COVID19_CaseStudy/data_acquisition
WORKSPACE_ROOT_GUESS="$( dirname "${SCRIPT_DIR}" )" # Should be DIVAS_COVID19_CaseStudy
ARRAYEXPRESS_DATA_DIR="${SCRIPT_DIR}/arrayexpress_data"
JSON_FILE="${ARRAYEXPRESS_DATA_DIR}/processed-data_filelist.json"
BASE_FTP_URL="ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/357/E-MTAB-9357/Files"
LOG_FILE="${ARRAYEXPRESS_DATA_DIR}/download_arrayexpress.log"

# Check if wget is installed
if ! command -v wget &> /dev/null
then
    echo "wget could not be found. Please install wget to use this script."
    exit 1
fi

# Check if python3 is installed (for parsing JSON)
if ! command -v python3 &> /dev/null
then
    echo "python3 could not be found. Please install python3 to use this script."
    exit 1
fi

# Check if the JSON file exists
if [ ! -f "${JSON_FILE}" ]; then
    echo "ERROR: JSON file list not found at ${JSON_FILE}"
    echo "Please ensure processed-data_filelist.json is downloaded to ${ARRAYEXPRESS_DATA_DIR}"
    exit 1
fi

echo "Starting download of ArrayExpress files for E-MTAB-9357..." | tee -a "${LOG_FILE}"
echo "Log file: ${LOG_FILE}" | tee -a "${LOG_FILE}"
echo "Target directory: ${ARRAYEXPRESS_DATA_DIR}" | tee -a "${LOG_FILE}"
echo "--------------------------------------------------" | tee -a "${LOG_FILE}"

# Create target directory if it doesn't exist
mkdir -p "${ARRAYEXPRESS_DATA_DIR}"

# Use Python to parse JSON and generate wget commands
python3 -c "
import json
import os

json_file_path = '${JSON_FILE}'
base_url = '${BASE_FTP_URL}'
target_dir = '${ARRAYEXPRESS_DATA_DIR}'
log_file = '${LOG_FILE}'

try:
    with open(json_file_path, 'r') as f:
        files_to_download = json.load(f)
except json.JSONDecodeError as e:
    print(f'Error decoding JSON: {e}')
    with open(log_file, 'a') as lf:
        lf.write(f'Error decoding JSON: {e}\n')
    exit(1)
except FileNotFoundError:
    print(f'Error: JSON file not found at {json_file_path}')
    with open(log_file, 'a') as lf:
        lf.write(f'Error: JSON file not found at {json_file_path}\n')
    exit(1)

downloaded_count = 0
skipped_count = 0
failed_count = 0
total_files = len(files_to_download)

print(f'Found {total_files} files to potentially download.')
with open(log_file, 'a') as lf:
    lf.write(f'Found {total_files} files to potentially download.\n')

for i, file_info in enumerate(files_to_download):
    file_path = file_info.get('path')
    if not file_path:
        print(f'Skipping entry with no path: {file_info}')
        with open(log_file, 'a') as lf:
            lf.write(f'Skipping entry with no path: {file_info}\n')
        continue

    file_name = os.path.basename(file_path)
    full_url = f'{base_url}/{file_path}'
    local_file_path = os.path.join(target_dir, file_name)
    
    # Use wget's timestamping and continue features.
    # -N: turn on timestamping. Don't re-retrieve files unless newer than local.
    # -c: continue getting a partially-downloaded file.
    # -P: save files to target_dir
    # -nv: non-verbose
    command = f'wget -N -c -P \'{target_dir}\' \'{full_url}\''
    
    print(f'Downloading ({i+1}/{total_files}): {file_name}...')
    with open(log_file, 'a') as lf:
        lf.write(f'Attempting to download: {full_url} to {local_file_path}\n')
        lf.write(f'Command: {command}\n')
        
    # Before running, check if file exists and has non-zero size to decide if it was already downloaded
    # This is a more robust check than relying purely on wget -N for skipping if file sizes match
    # However, for simplicity with wget -N and -c, we let wget handle this.
    
    exit_code = os.system(command)
    
    if exit_code == 0:
        print(f'Successfully downloaded or updated {file_name}')
        with open(log_file, 'a') as lf:
            lf.write(f'Successfully downloaded/updated {file_name}\n')
        downloaded_count +=1 # Technically, this counts success, which might be an update or skip due to -N
    else:
        print(f'Failed to download {file_name}. Exit code: {exit_code}')
        with open(log_file, 'a') as lf:
            lf.write(f'Failed to download {file_name}. Exit code: {exit_code}\n')
        failed_count += 1

print(f'\nDownload process completed.')
print(f'Total files processed: {total_files}')
print(f'Successfully downloaded/updated: {downloaded_count}')
print(f'Failed: {failed_count}')

with open(log_file, 'a') as lf:
    lf.write('\nDownload process completed.\n')
    lf.write(f'Total files processed: {total_files}\n')
    lf.write(f'Successfully downloaded/updated: {downloaded_count}\n')
    lf.write(f'Failed: {failed_count}\n')

" | tee -a "${LOG_FILE}"

echo "--------------------------------------------------" | tee -a "${LOG_FILE}"
echo "ArrayExpress file download script finished. Check log: ${LOG_FILE}" | tee -a "${LOG_FILE}" 