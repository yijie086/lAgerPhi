#!/bin/bash
# This script has two modes:
#
# 1) Launcher mode (no arguments):
#    - Splits the full run range into chunks (<= MAX_CHUNK runs each)
#    - Submits each chunk as a background job using nohup
#
#       ./run_lund_large.sh
#
# 2) Worker mode (two arguments: START_RUN END_RUN):
#    - Actually runs bin/lager, calls clean_lund, and deletes temporary files
#      for runs in [START_RUN, END_RUN].
#
#       ./run_lund_large.sh 1 100
#
# All configuration (NEV, TOTAL_RUNS, etc.) is set below.
#
# This version:
#   - Can be launched from ANY directory
#   - Keeps all outputs in a separate output directory
#   - Still uses launcher/worker mode with nohup

# ---------------- Paths & configuration ----------------

# Directory where THIS script lives (i.e. your lAgerPhi dir)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_NAME="$(basename "$0")"
SCRIPT_PATH="${SCRIPT_DIR}/${SCRIPT_NAME}"

# Root of lAger installation (by default, same directory as this script)
LAGER_HOME="${LAGER_HOME:-${SCRIPT_DIR}}"

# JSON configuration file for lAgerPhi (adjust name if needed)
CONFIG="${LAGER_HOME}/CLAS1210GeV.ep-phi.gen.json"

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

# Output directory:
# - Default: ./lager_runs relative to WHERE YOU LAUNCH the script
# - Can override by: OUT_DIR=/some/path ./run_lund_large.sh
OUT_DIR="${OUT_DIR:-${PWD}/lager_runs}"
mkdir -p "${OUT_DIR}"

# ---------------- Launcher mode ----------------
# No arguments → submit worker jobs with nohup
if [ "$#" -eq 0 ]; then
    echo "Launcher mode: submitting nohup jobs."
    echo "  Total runs   : ${TOTAL_RUNS}"
    echo "  Chunk size   : ${MAX_CHUNK}"
    echo "  Events / run : ${NEV}"
    echo "  Output dir   : ${OUT_DIR}"
    echo

    start=1
    while [ "$start" -le "$TOTAL_RUNS" ]; do
        end=$(( start + MAX_CHUNK - 1 ))
        if [ "$end" -gt "$TOTAL_RUNS" ]; then
            end="$TOTAL_RUNS"
        fi

        LOG_FILE="lager_${start}_${end}.log"
        echo "Submitting runs ${start} to ${end} -> log: ${LOG_FILE}"

        # Export OUT_DIR so worker sees the same directory
        OUT_DIR="${OUT_DIR}" nohup "${SCRIPT_PATH}" "$start" "$end" > "${LOG_FILE}" 2>&1 &

        start=$(( end + 1 ))
    done

    echo
    echo "All nohup jobs submitted."
    exit 0
fi

# ---------------- Worker mode ----------------
# Two arguments: START_RUN END_RUN
if [ "$#" -ne 2 ]; then
    echo "Usage:"
    echo "  Launcher mode: ./${SCRIPT_NAME}"
    echo "  Worker mode  : ./${SCRIPT_NAME} START_RUN END_RUN"
    exit 1
fi

START_RUN="$1"
END_RUN="$2"

echo "Worker mode: processing runs ${START_RUN} to ${END_RUN}"
echo "CONFIG   = ${CONFIG}"
echo "NEV      = ${NEV}"
echo "PREFIX   = ${PREFIX}"
echo "OUT_DIR  = ${OUT_DIR}"
echo "LAGER_HOME = ${LAGER_HOME}"
echo

for ((i=START_RUN; i<=END_RUN; i++)); do
    # Format run number as 00001, 00002, ...
    RUN_TAG=$(printf "%05d" "$i")

    echo "=== Running lager for run ${RUN_TAG} ==="

    # Run generator, output to OUT_DIR
    "${LAGER_HOME}/bin/lager" -c "${CONFIG}" -r "${i}" -e "${NEV}" -o "${OUT_DIR}"

    # Base name for this run, without extension
    # Example: CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000
    BASE="${PREFIX}.run${RUN_TAG}-${NEV}"

    # Expected .gemc LUND file (inside OUT_DIR)
    IN_LUND="${OUT_DIR}/${BASE}.gemc"

    # Cleaned output file (txt) (also in OUT_DIR)
    OUT_TXT="${OUT_DIR}/${BASE}_clean.txt"

    # Clean if LUND exists
    if [ -f "${IN_LUND}" ]; then
        echo "Cleaning ${IN_LUND} -> ${OUT_TXT}"
        "${LAGER_HOME}/clean_lund" "${IN_LUND}" "${OUT_TXT}"
    else
        echo "WARNING: LUND file not found: ${IN_LUND}"
    fi

    echo "Deleting temporary files for run ${RUN_TAG} ..."

    # Delete .gemc, .log, .root, .json for this run (from OUT_DIR)
    for ext in gemc log root json; do
        FILE="${OUT_DIR}/${BASE}.${ext}"
        if [ -f "${FILE}" ]; then
            rm -f "${FILE}"
            echo "  Deleted: ${FILE}"
        fi
    done

    echo
done

echo "Worker finished: runs ${START_RUN} to ${END_RUN}."
