#!/bin/bash

# Script to organize downloaded ArrayExpress files from E-MTAB-9357 into subdirectories based on omics type.

# Ensure the script is run from the workspace root or adjust paths accordingly
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # Should be DIVAS-code/DATA_Covid/downloads
WORKSPACE_ROOT_GUESS="$( dirname "$( dirname "${SCRIPT_DIR}" )" )" # Should be DIVAS-code
ARRAYEXPRESS_BASE_DIR="${WORKSPACE_ROOT_GUESS}/DATA_Covid/downloads/arrayexpress_data"

# Define target directories
GEX_DIR="${ARRAYEXPRESS_BASE_DIR}/gex_data"
PRO_DIR="${ARRAYEXPRESS_BASE_DIR}/pro_data"
CD8_TCR_DIR="${ARRAYEXPRESS_BASE_DIR}/cd8_tcr_data"
CD4_TCR_DIR="${ARRAYEXPRESS_BASE_DIR}/cd4_tcr_data"
BCR_DIR="${ARRAYEXPRESS_BASE_DIR}/bcr_data"

LOG_FILE="${ARRAYEXPRESS_BASE_DIR}/organize_files.log"

echo "Starting organization of ArrayExpress files..." | tee -a "${LOG_FILE}"
echo "Source directory: ${ARRAYEXPRESS_BASE_DIR}" | tee -a "${LOG_FILE}"
echo "Log file: ${LOG_FILE}" | tee -a "${LOG_FILE}"
echo "--------------------------------------------------" | tee -a "${LOG_FILE}"

# Ensure target directories exist (should have been created by a previous step)
mkdir -p "${GEX_DIR}" "${PRO_DIR}" "${CD8_TCR_DIR}" "${CD4_TCR_DIR}" "${BCR_DIR}"

shopt -s nullglob # Ensures that if no files match a pattern, the loop doesn't run with the pattern itself

# Counters
moved_gex=0
moved_pro=0
moved_cd8_tcr=0
moved_cd4_tcr=0
moved_bcr=0
unmatched_files=0

cd "${ARRAYEXPRESS_BASE_DIR}" || exit

for file in *.txt.gz; do
    if [[ "$file" == *"_gex_library_"* ]]; then
        mv -- "$file" "${GEX_DIR}/" && echo "Moved to gex_data: $file" | tee -a "${LOG_FILE}"
        ((moved_gex++))
    elif [[ "$file" == *"_pro_library_"* ]]; then
        mv -- "$file" "${PRO_DIR}/" && echo "Moved to pro_data: $file" | tee -a "${LOG_FILE}"
        ((moved_pro++))
    elif [[ "$file" == *"_cd8_tcr_library_"* ]]; then
        mv -- "$file" "${CD8_TCR_DIR}/" && echo "Moved to cd8_tcr_data: $file" | tee -a "${LOG_FILE}"
        ((moved_cd8_tcr++))
    elif [[ "$file" == *"_cd4_tcr_library_"* ]]; then
        mv -- "$file" "${CD4_TCR_DIR}/" && echo "Moved to cd4_tcr_data: $file" | tee -a "${LOG_FILE}"
        ((moved_cd4_tcr++))
    elif [[ "$file" == *"_bcr_library_"* ]]; then
        mv -- "$file" "${BCR_DIR}/" && echo "Moved to bcr_data: $file" | tee -a "${LOG_FILE}"
        ((moved_bcr++))
    else
        echo "WARNING: Unmatched file: $file" | tee -a "${LOG_FILE}"
        ((unmatched_files++))
    fi
done

cd "${WORKSPACE_ROOT_GUESS}" || exit # Go back to original directory if needed

echo "--------------------------------------------------" | tee -a "${LOG_FILE}"
echo "File organization completed." | tee -a "${LOG_FILE}"
echo "Summary:" | tee -a "${LOG_FILE}"
echo "Moved to gex_data: ${moved_gex}" | tee -a "${LOG_FILE}"
echo "Moved to pro_data: ${moved_pro}" | tee -a "${LOG_FILE}"
echo "Moved to cd8_tcr_data: ${moved_cd8_tcr}" | tee -a "${LOG_FILE}"
echo "Moved to cd4_tcr_data: ${moved_cd4_tcr}" | tee -a "${LOG_FILE}"
echo "Moved to bcr_data: ${moved_bcr}" | tee -a "${LOG_FILE}"
echo "Unmatched files: ${unmatched_files}" | tee -a "${LOG_FILE}"

TOTAL_MOVED=$((moved_gex + moved_pro + moved_cd8_tcr + moved_cd4_tcr + moved_bcr))
echo "Total files moved: ${TOTAL_MOVED}" | tee -a "${LOG_FILE}"

if [ "${unmatched_files}" -ne 0 ]; then
    echo "WARNING: ${unmatched_files} file(s) could not be categorized. Please check the log and the source directory." | tee -a "${LOG_FILE}"
fi

if [ "${TOTAL_MOVED}" -eq 1340 ] && [ "${unmatched_files}" -eq 0 ]; then
    echo "All 1340 files successfully categorized and moved." | tee -a "${LOG_FILE}"
else
    echo "Please review the log for details on categorization and any unmatched files." | tee -a "${LOG_FILE}"
fi 