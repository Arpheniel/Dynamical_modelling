#!/bin/bash
# Parallel grid launcher for PGC 35706 around best-fit (RHALO=107, VHALO=154).
#
# BEFORE RUNNING: activate the agama virtualenv!
#     agamaenv
#   (which expands to: source ~/.environments/3_11/bin/activate)
#
# Stage 1: 5x5 grid at N=10000 (25 points, one already exists, so 24 new).
# Stage 2: 3 cross-checks at N=40000.
#
# Launches up to $PARALLEL points concurrently. Skips points whose .npz
# already exists. Results are appended to resultsGH.txt; each per-point
# log goes to logs/<stem>.log so parallel processes do not interleave
# their output.
#
# Usage:
#   agamaenv
#   ./run_grid_PGC.sh              # defaults to PARALLEL=6
#   PARALLEL=10 ./run_grid_PGC.sh  # override concurrency
#
# Monitoring from another ssh session:
#   tail -f grid_log.txt                   # high-level progress
#   tail -f logs/M1e+07*Rh117_Vh154*.log   # a specific point
#   ls -1 *_N10000_*.npz | wc -l           # how many points are done
#
# To stop early: kill parent bash (Ctrl+C) then
#   pkill -f run_forstand_grid_PGC_35706.py

set -u
PY=python3
SCRIPT=run_forstand_grid_PGC_35706.py
PARALLEL=${PARALLEL:-6}
LOGDIR=logs
MASTER_LOG=grid_log.txt
mkdir -p "$LOGDIR"

# sanity check: are we in the right virtualenv?
if ! python3 -c "import agama" >/dev/null 2>&1; then
    echo "ERROR: 'import agama' failed.  Did you run 'agamaenv' first?" >&2
    exit 1
fi

echo "====== Grid run started at $(date) ======" | tee -a "$MASTER_LOG"
echo "Parallelism: $PARALLEL concurrent processes" | tee -a "$MASTER_LOG"

# Throttle: wait until fewer than $PARALLEL of our jobs are running.
throttle () {
    while (( $(jobs -rp | wc -l) >= PARALLEL )); do
        wait -n 2>/dev/null || sleep 2
    done
}

submit () {
    local RH=$1 VH=$2 N=$3 SAVE=$4
    local STEM="M1e+07_O0_Rh${RH}_Vh${VH}_i42_a0_N${N}_R0.00_GH_DensitySphHarm"
    if ls "${STEM}".npz >/dev/null 2>&1; then
        echo "[skip] ${STEM}.npz exists" | tee -a "$MASTER_LOG"
        return 0
    fi
    throttle
    local LOG="${LOGDIR}/${STEM}.log"
    (
        echo "[run ] $STEM  start $(date +%H:%M:%S)" >> "$MASTER_LOG"
        $PY "$SCRIPT" RHALO=$RH VHALO=$VH NUMORBITS=$N SAVE_ORBITS=$SAVE \
            > "$LOG" 2>&1
        rc=$?
        if [[ $rc -eq 0 ]]; then
            echo "[done] $STEM  end   $(date +%H:%M:%S)" >> "$MASTER_LOG"
        else
            echo "[FAIL] $STEM  rc=$rc (see $LOG)" >> "$MASTER_LOG"
        fi
    ) &
}

###############################################################
# Stage 1: 5x5 grid at N=10000 around (107, 154)
###############################################################
echo "--- Stage 1: N=10000 grid (up to 24 new points) ---" | tee -a "$MASTER_LOG"
for RH in 87 97 107 117 127; do
    for VH in 124 139 154 169 184; do
        submit $RH $VH 10000 no
    done
done
wait
echo "--- Stage 1 complete at $(date) ---" | tee -a "$MASTER_LOG"

###############################################################
# Stage 2: N=40000 cross-checks (heavier per-point, ~4x)
###############################################################
echo "--- Stage 2: N=40000 stability check ---" | tee -a "$MASTER_LOG"
submit 117 154 40000 no
submit 107 169 40000 no
submit 127 169 40000 no
wait

echo "====== Grid run finished at $(date) ======" | tee -a "$MASTER_LOG"
echo "resultsGH.txt lines: $(wc -l < resultsGH.txt 2>/dev/null || echo 0)" | tee -a "$MASTER_LOG"
