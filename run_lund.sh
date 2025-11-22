#!/bin/bash
# This script automatically:
#   1. Runs bin/lager for multiple runs
#   2. Cleans the generated LUND files using clean_lund
#   3. Deletes .gemc, .log, and .root files after cleaning
#
# Modify NRUNS / NEV / CONFIG / PREFIX as needed.

# JSON configuration file for lAgerPhi
CONFIG="CLAS1210GeV.ep-phi.gen.json"

# Number of events per run (defines the final "-1000" part)
NEV=5000

# Total number of runs to generate
NRUNS=10

# Prefix of the lAger output files
# Example: CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-1000.gemc
PREFIX="CLAS-ep-phi-10GeV.ep-phi.4pi"

for ((i=1; i<=NRUNS; i++)); do
    # Format run number as 00001, 00002, ...
    RUN_TAG=$(printf "%05d" "$i")

    echo "=== Running lager for run ${RUN_TAG} ==="

    # Run generator
    bin/lager -c "${CONFIG}" -r "${i}" -e "${NEV}" -o .

    # Expected .gemc LUND file
    IN_LUND="${PREFIX}.run${RUN_TAG}-${NEV}.gemc"

    # Cleaned output file (txt)
    OUT_TXT="${PREFIX}.run${RUN_TAG}-${NEV}_clean.txt"

    # Clean if exists
    if [ -f "${IN_LUND}" ]; then
        echo "Cleaning ${IN_LUND} -> ${OUT_TXT}"
        ./clean_lund "${IN_LUND}" "${OUT_TXT}"
    else
        echo "WARNING: LUND file not found: ${IN_LUND}"
    fi

    echo "Deleting temporary files..."

    # Delete .gemc file
    if [ -f "${IN_LUND}" ]; then
        rm -f "${IN_LUND}"
        echo "  Deleted: ${IN_LUND}"
    fi

    # Delete .log (produced by lager)
    LOG_FILE="${PREFIX}.run${RUN_TAG}-${NEV}.log"
    if [ -f "${LOG_FILE}" ]; then
        rm -f "${LOG_FILE}"
        echo "  Deleted: ${LOG_FILE}"
    fi

    # Delete any .root file created during the run
    ROOT_FILE="${PREFIX}.run${RUN_TAG}-${NEV}.root"
    if [ -f "${ROOT_FILE}" ]; then
        rm -f "${ROOT_FILE}"
        echo "  Deleted: ${ROOT_FILE}"
    fi

    JSON_FILE="${PREFIX}.run${RUN_TAG}-${NEV}.json"
    if [ -f "${JSON_FILE}" ]; then
        rm -f "${JSON_FILE}"
        echo "  Deleted: ${JSON_FILE}"
    fi

    LOG_FILE="${PREFIX}.run${RUN_TAG}.log"
    if [ -f "${LOG_FILE}" ]; then
        rm -f "${LOG_FILE}"
        echo "  Deleted: ${LOG_FILE}"
    fi

    echo
done

echo "All runs completed."
