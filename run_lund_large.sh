#!/bin/bash
# run_lund_large.sh (robust version)

set -euo pipefail

# ---------------- Paths & configuration ----------------

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_NAME="$(basename "$0")"
SCRIPT_PATH="${SCRIPT_DIR}/${SCRIPT_NAME}"

LAGER_HOME="${LAGER_HOME:-${SCRIPT_DIR}}"
CONFIG="${LAGER_HOME}/CLAS1210GeV.ep-phi.gen.json"
PREFIX="CLAS-ep-phi-10GeV.ep-phi.4pi"

NEV=5000
TOTAL_RUNS=10
MAX_CHUNK=1

FMAX=1.0
TARGET_NUM=9
TARGET_DEN=10

OUT_DIR="${OUT_DIR:-${PWD}/lager_runs}"
mkdir -p "${OUT_DIR}"

SUMMARY="${OUT_DIR}/sampling_summary.txt"

usage() {
  echo "Usage:"
  echo "  Launcher mode: ./${SCRIPT_NAME}"
  echo "  Worker mode  : ./${SCRIPT_NAME} START_RUN END_RUN"
}

# ---------------- Launcher mode ----------------
if [ "$#" -eq 0 ]; then
  echo "Launcher mode: submitting nohup jobs."
  echo "  Total runs   : ${TOTAL_RUNS}"
  echo "  Chunk size   : ${MAX_CHUNK}"
  echo "  Events/batch : ${NEV}"
  echo "  Target       : ${TARGET_NUM}/${TARGET_DEN} of NEV"
  echo "  Output dir   : ${OUT_DIR}"
  echo "  Summary file : ${SUMMARY}"
  echo

  : > "${SUMMARY}"

  start=1
  while [ "${start}" -le "${TOTAL_RUNS}" ]; do
    end=$(( start + MAX_CHUNK - 1 ))
    if [ "${end}" -gt "${TOTAL_RUNS}" ]; then
      end="${TOTAL_RUNS}"
    fi

    LOG_FILE="${OUT_DIR}/lager_${start}_${end}.log"
    echo "Submitting runs ${start} to ${end} -> log: ${LOG_FILE}"

    # Force bash explicitly (avoids any shebang / nohup weirdness)
    OUT_DIR="${OUT_DIR}" \
    LAGER_HOME="${LAGER_HOME}" \
    nohup /bin/bash "${SCRIPT_PATH}" "${start}" "${end}" > "${LOG_FILE}" 2>&1 &

    start=$(( end + 1 ))
  done

  echo
  echo "All nohup jobs submitted."
  exit 0
fi

# ---------------- Worker mode ----------------
if [ "$#" -ne 2 ]; then
  usage
  exit 1
fi

# Debug trace in worker logs
# set -x

START_RUN="$1"
END_RUN="$2"

echo "Worker mode: processing runs ${START_RUN} to ${END_RUN}"
echo "CONFIG      = ${CONFIG}"
echo "NEV         = ${NEV}"
echo "PREFIX      = ${PREFIX}"
echo "OUT_DIR     = ${OUT_DIR}"
echo "LAGER_HOME  = ${LAGER_HOME}"
echo "FMAX        = ${FMAX}"
echo "TARGET      = ${TARGET_NUM}/${TARGET_DEN}"
echo "SUMMARY     = ${SUMMARY}"
echo

touch "${SUMMARY}"

run="${START_RUN}"
while [ "${run}" -le "${END_RUN}" ]; do
  RUN_TAG="$(printf "%05d" "${run}")"
  FINAL="${OUT_DIR}/run${RUN_TAG}_final.txt"
  : > "${FINAL}"

  TARGET_ACC=$(( (NEV * TARGET_NUM) / TARGET_DEN ))

  ACC=0
  BATCH=0

  echo "=== RUN ${RUN_TAG} start (target accepted >= ${TARGET_ACC}) ==="

  while [ "${ACC}" -lt "${TARGET_ACC}" ]; do
    BATCH=$((BATCH + 1))

    # Unique deterministic seed per batch
    SEED=$(( run*1000000 + BATCH ))
    SEED_TAG="$(printf "%05d" "${SEED}")"

    echo "Run ${RUN_TAG} batch ${BATCH} (seed=${SEED_TAG})"

    # Generate
    "${LAGER_HOME}/bin/lager" -c "${CONFIG}" -r "${SEED}" -e "${NEV}" -o "${OUT_DIR}" -v 0

    # lager output base uses seed (because -r embedded in filename)
    BASE="${PREFIX}.run${SEED_TAG}-${NEV}"
    IN_LUND="${OUT_DIR}/${BASE}.gemc"
    OUT_TMP="${OUT_DIR}/tmp_run${RUN_TAG}_batch${BATCH}.lund"

    if [ ! -f "${IN_LUND}" ]; then
      echo "[ERROR] Missing LUND file: ${IN_LUND}"
      exit 3
    fi

    # Clean + rejection sampling (your interface)
    RES=$("${LAGER_HOME}/clean_lund" \
      "${IN_LUND}" "${OUT_TMP}" \
      "${FMAX}" "${SUMMARY}" \
      "run${RUN_TAG}_batch${BATCH}")

    ADD="$(echo "${RES}" | sed 's/ACCEPTED_EVENTS=//')"
    case "${ADD}" in
      (*[!0-9]*|'') ADD=0 ;;
    esac

    ACC=$((ACC + ADD))
    echo "  Accepted this batch: ${ADD} | Total accepted: ${ACC}"

    cat "${OUT_TMP}" >> "${FINAL}"

    # cleanup
    rm -f "${IN_LUND}" "${OUT_TMP}"
    for ext in log root json; do
      f="${OUT_DIR}/${BASE}.${ext}"
      [ -f "${f}" ] && rm -f "${f}"
    done
  done

  echo "=== RUN ${RUN_TAG} done: accepted ${ACC} (saved: ${FINAL}) ==="
  echo

  run=$((run + 1))
done

echo "Worker finished: runs ${START_RUN} to ${END_RUN}."
