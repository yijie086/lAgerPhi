#!/bin/bash

CONFIG="CLAS1210GeV.ep-phi.gen.json"
PREFIX="CLAS-ep-phi-10GeV.ep-phi.4pi"

NEV_GEN=10000
TARGET=$NEV_GEN          # <-- 修正：原来写成 TARGET=NEV_GEN 会变成字符串
NRUNS=3
FMAX=1.0

SUMMARY="sampling_summary.txt"
> "$SUMMARY"

for ((r=1; r<=NRUNS; r++)); do

  RUN=$(printf "%05d" "$r")
  FINAL="run${RUN}_final.txt"
  > "$FINAL"

  ACC=0
  BATCH=0

  while (( ACC < TARGET*9/10 )); do
    ((BATCH++))

    # -----------------------------
    # Deterministic, reproducible seed per batch
    # -----------------------------
    SEED=$(( r*1000000 + BATCH ))
    SEED_TAG=$(printf "%05d" "$SEED")

    echo "Run $RUN batch $BATCH (seed=$SEED_TAG)"

    # Generate with unique seed each batch
    bin/lager -c "$CONFIG" -r "$SEED" -e "$NEV_GEN" -o . -v 0

    # IMPORTANT: lager output file name follows the -r value (seed)
    IN="${PREFIX}.run${SEED_TAG}-${NEV_GEN}.gemc"
    OUT="tmp_${RUN}_${BATCH}.lund"

    # Clean + rejection sampling
    RES=$(./clean_lund \
      "$IN" "$OUT" \
      "$FMAX" "$SUMMARY" \
      "run${RUN}_batch${BATCH}")

    ADD=$(echo "$RES" | sed 's/ACCEPTED_EVENTS=//')
    ACC=$((ACC + ADD))

    cat "$OUT" >> "$FINAL"

    # Clean intermediate files for this batch/seed
    rm -f "$IN" "$OUT"

    # Delete .log/.root/.json produced by lager for this seed
    for ext in log root json; do
      f="${PREFIX}.run${SEED_TAG}-${NEV_GEN}.${ext}"
      if [ -f "$f" ]; then
        rm -f "$f"
        echo "  Deleted: $f"
      fi
    done
  done

  echo "Run $RUN done, accepted $ACC"

done
