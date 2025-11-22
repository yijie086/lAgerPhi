#!/bin/bash
# This script has two modes:
#
# 1) Launcher mode (no arguments):
#    - Splits the full run range into chunks (<= MAX_CHUNK runs each)
#    - Submits each chunk as a background job using nohup
#
#       ./run_lager_clean.sh
#
# 2) Worker mode (two arguments: START_RUN END_RUN):
#    - Actually runs bin/lager, calls clean_lund, and deletes temporary files
#      for runs in [START_RUN, END_RUN].
#
#       ./run_lager_clean.sh 1 100
#
# All configuration (NEV, TOTAL_RUNS, etc.) is set below.

# ---------------- General configuration ----------------

# JSON configuration file for lAgerPhi
CONFIG="CLAS1210GeV.ep-phi.gen.json"

# Number of events per run (defines the final "-5000" part in filenames)
NEV=5000

# Total number of runs you want to generate in total (for launcher mode)
TOTAL_RUNS=100

# Maximum number of runs per background job
MAX_CHUNK=10

# Prefix of the lAger output files
# Example file name pattern:
#   CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000.gemc
PREFIX="CLAS-ep-phi-10GeV.ep-phi.4pi"

# Name of this script (used for recursive calls)
SCRIPT_NAME="$(basename "$0")"

# ------------------------------------------------------
# Launcher mode: no or one argument (we only use TOTAL_RUNS)
# ------------------------------------------------------
if [ "$#" -eq 0 ]; then
    echo "Launcher mode: submitting nohup jobs."
    echo "  Total runs   : ${TOTAL_RUNS}"
    echo "  Chunk size   : ${MAX_CHUNK}"
    echo "  Events / run : ${NEV}"
    echo

    start=1
    while [ "$start" -le "$TOTAL_RUNS" ]; do
        end=$(( start + MAX_CHUNK - 1 ))
        if [ "$end" -gt "$TOTAL_RUNS" ]; then
            end="$TOTAL_RUNS"
        fi

        LOG_FILE="lager_${start}_${end}.log"
        echo "Submitting runs ${start} to ${end} -> log: ${LOG_FILE}"

        # Submit this script in worker mode as a background job via nohup
        nohup "./${SCRIPT_NAME}" "$start" "$end" > "${LOG_FILE}" 2>&1 &

        start=$(( end + 1 ))
    done

    echo
    echo "All nohup jobs submitted."
    exit 0
fi

# ------------------------------------------------------
# Worker mode: two arguments (START_RUN END_RUN)
# ------------------------------------------------------
if [ "$#" -ne 2 ]; then
    echo "Usage:"
    echo "  Launcher mode: ./${SCRIPT_NAME}"
    echo "  Worker mode  : ./${SCRIPT_NAME} START_RUN END_RUN"
    exit 1
fi

START_RUN="$1"
END_RUN="$2"

echo "Worker mode: processing runs ${START_RUN} to ${END_RUN}"
echo "CONFIG = ${CONFIG}, NEV = ${NEV}, PREFIX = ${PREFIX}"
echo

for ((i=START_RUN; i<=END_RUN; i++)); do
    # Format run number as 00001, 00002, ...
    RUN_TAG=$(printf "%05d" "$i")

    echo "=== Running lager for run ${RUN_TAG} ==="

    # Run generator, output to current directory
    bin/lager -c "${CONFIG}" -r "${i}" -e "${NEV}" -o .

    # Base name for this run, without extension
    # Example: CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000
    BASE="${PREFIX}.run${RUN_TAG}-${NEV}"

    # Expected .gemc LUND file
    IN_LUND="${BASE}.gemc"

    # Cleaned output file (txt)
    OUT_TXT="${BASE}_clean.txt"

    # Clean if LUND exists
    if [ -f "${IN_LUND}" ]; then
        echo "Cleaning ${IN_LUND} -> ${OUT_TXT}"
        ./clean_lund "${IN_LUND}" "${OUT_TXT}"
    else
        echo "WARNING: LUND file not found: ${IN_LUND}"
    fi

    echo "Deleting temporary files for run ${RUN_TAG} ..."

    # Delete .gemc, .log, .root, .json for this run
    for ext in gemc log root json; do
        FILE="${BASE}.${ext}"
        if [ -f "${FILE}" ]; then
            rm -f "${FILE}"
            echo "  Deleted: ${FILE}"
        fi
    done

    echo
done

echo "Worker finished: runs ${START_RUN} to ${END_RUN}."
