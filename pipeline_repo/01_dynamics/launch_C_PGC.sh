#!/bin/bash
# launch_C_PGC.sh -- APPROACH C (1-D mass-normalization scan at a
# cosmologically fixed NFW scale radius) for PGC 35706.
# Place next to forstand_C_PGC_35706.py and the galaxy's data files.
# Run:  bash launch_C_PGC.sh
#
# Approach C fixes r_s at the abundance-matching value (computed inside
# the python script from Moster+13 + Dutton&Maccio14) and scans only the
# halo normalization v_halo (NVHALO=12 points by default), with Upsilon
# marginalized per point by the sqrt(Y) trick. So the whole job is
#   12 orbital libraries (one QP-solve each) + 1 final orbit-saving solve.
# It is fully disjoint from the _R1 2-D grid (separate _C1 output names),
# so it can run in the same folder without touching the R1 files.
#
# Workflow:
#   STEP 1 -- scan:      4 parallel workers (chunks) share the 12 v_halo points
#   STEP 2 -- postgrid:  merge + verify + final orbit-saving run, 1 process
#
# Outputs (all _C1-tagged):
#   grid_chi2_GH_C1.npz           -- master chi^2(v_halo) curve (1 x NVHALO)
#   grid_chi2_GH_C1_chunk*.npz    -- per-chunk fragments (merged at postgrid)
#   resultsGH_C1.txt              -- per-Y line-summary from the final run
#   resultsGH_C1_chunk*.txt       -- per-chunk per-Y line-summary
#   resultsGH_C1_verify.txt       -- verify (Y baked, global-rescale)
#   verify_GH_C1.npz              -- verify table around the scan minimum
#   M*_..._GH_*.nbody             -- final orbit-resolved model (the deliverable)
#   logs/chunkC*_PGC.log          -- per-worker stdout/stderr
#   logs/postgridC_PGC.log        -- postgrid stdout/stderr
#
# Auto-resume:
#   If killed mid-scan, just rerun `bash launch_C_PGC.sh`. Per-cell
#   subprocess isolation rewrites grid_chi2_GH_C1.npz after every
#   successful point, so the rerun skips already-done points.

set -e
set -o pipefail

# ---- settings ----
NCHUNKS=4
THREADS_PER_CHUNK="${THREADS_PER_CHUNK:-5}"        # 4 x 5 = 20 threads while scanning
THREADS_TOTAL="${THREADS_TOTAL:-20}"               # postgrid single-process threads

# ---- pick the python script by folder content ----
PY=$(ls forstand_C_PGC_*.py 2>/dev/null | head -1)
if [ -z "$PY" ]; then
    echo "ERROR: no forstand_C_PGC_*.py in current folder ($(pwd))"
    exit 1
fi

mkdir -p logs

echo "============================================================"
echo "  Approach C (1-D v_halo scan, fixed cosmological r_s)"
echo "  Galaxy folder        : $(pwd)"
echo "  Python script        : $PY"
echo "  Scan: $NCHUNKS workers x $THREADS_PER_CHUNK threads"
echo "  Postgrid: 1 process  x $THREADS_TOTAL threads"
echo "============================================================"

# ============================================================
#  STEP 1 -- parallel 1-D scan
# ============================================================
echo ""
echo ">>> STEP 1: parallel v_halo scan (C)"
echo ""

export OMP_NUM_THREADS=$THREADS_PER_CHUNK
export PYTHONFAULTHANDLER=1

declare -a PIDS
for i in $(seq 0 $((NCHUNKS-1))); do
    LOG="logs/chunkC${i}_PGC.log"
    python -u "$PY" phase=grid chunk=$i nchunks=$NCHUNKS > "$LOG" 2>&1 &
    PIDS[$i]=$!
    echo "  chunk $i -> PID=${PIDS[$i]}  log=$LOG"
done

echo ""
echo "All $NCHUNKS workers launched. OMP_NUM_THREADS=$THREADS_PER_CHUNK each."
echo "(detach screen with Ctrl-A D; the workers keep running.)"
echo ""

# ---- wait for each PID individually and check exit codes ----
set +e
N_FAILED=0
declare -a FAILED_CHUNKS
for i in $(seq 0 $((NCHUNKS-1))); do
    pid=${PIDS[$i]}
    if wait "$pid"; then
        echo "  chunk $i (PID $pid): exit 0  OK"
    else
        rc=$?
        echo "  chunk $i (PID $pid): exit $rc  FAILED"
        N_FAILED=$((N_FAILED + 1))
        FAILED_CHUNKS+=("$i")
    fi
done
set -e

if [ $N_FAILED -gt 0 ]; then
    echo ""
    echo "============================================================"
    echo "  STEP 1 ABORTED: $N_FAILED / $NCHUNKS chunks failed: ${FAILED_CHUNKS[*]}"
    echo "  See logs/chunkC*_PGC.log for details."
    echo "  Re-run 'bash launch_C_PGC.sh' to auto-resume."
    echo "============================================================"
    exit 1
fi

echo ""
echo ">>> STEP 1 DONE: all $NCHUNKS workers finished cleanly"

# ============================================================
#  STEP 2 -- postgrid: merge + verify + final
# ============================================================
echo ""
echo ">>> STEP 2: postgrid (merge + verify + final)"
echo ""

export OMP_NUM_THREADS=$THREADS_TOTAL

python -u "$PY" phase=postgrid nchunks=$NCHUNKS 2>&1 | tee logs/postgridC_PGC.log

echo ""
echo "============================================================"
echo "  Approach C ALL DONE in $(pwd)"
echo "  Results:"
echo "    grid_chi2_GH_C1.npz       -- chi^2(v_halo) curve"
echo "    verify_GH_C1.npz          -- Y-baked global-rescale verification"
echo "    M*_N*_*.nbody             -- final orbit-resolved model"
echo "============================================================"
