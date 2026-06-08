#!/bin/bash
# finish_C_LEDA.sh -- STUBBORN finisher for the approach-C LEDA scan.
#
# Why this exists
# ---------------
# Two independent failure modes were biting the plain launch_C_LEDA.sh:
#   (1) the shared server randomly kills python -> a cell stays NaN;
#   (2) a Python-level exception inside a cell is SWALLOWED (onecell writes
#       NaN and exits 0), so the parent treats it as "done" and never retries.
# Either way a hole survives, and postgrid would then finalize on an
# incomplete grid (nanargmin silently skips the hole).
#
# This script loops the scan until the master grid is genuinely 12/12, then
# runs postgrid (which now also REFUSES to finalize on a holey grid and
# finalizes at the cosmological prior, FINAL_BY=prior).
#
# Usage:
#   bash finish_C_LEDA.sh                 # default: up to 30 scan passes
#   MAX_PASSES=60 bash finish_C_LEDA.sh   # more stubborn
#   FINAL_BY=chi2min bash finish_C_LEDA.sh# finalize at chi2 min instead of prior
# Safe to Ctrl-C and rerun: every pass auto-resumes (done cells are skipped).

set -u
set -o pipefail

NCHUNKS=4
THREADS_PER_CHUNK="${THREADS_PER_CHUNK:-5}"
THREADS_TOTAL="${THREADS_TOTAL:-20}"
MAX_PASSES="${MAX_PASSES:-30}"
SLEEP_BETWEEN="${SLEEP_BETWEEN:-20}"      # seconds to let the server settle between passes

PY=$(ls forstand_C_LEDA_*.py 2>/dev/null | head -1)
if [ -z "$PY" ]; then
    echo "ERROR: no forstand_C_LEDA_*.py in $(pwd)"; exit 1
fi
mkdir -p logs

export PYTHONFAULTHANDLER=1

# --- helper: how many master-grid cells are finite (prints "<finite> <total>") ---
grid_status() {
    python - <<'PYEOF' 2>/dev/null
import numpy, os
f='grid_chi2_GH_C1.npz'
if not os.path.exists(f):
    print("0 0"); raise SystemExit
import glob
# merge master + all chunk fragments the way postgrid will, just to count
chi=None
for g in sorted(set(['grid_chi2_GH_C1.npz']+glob.glob('grid_chi2_GH_C1_chunk*.npz'))):
    try: d=numpy.load(g)
    except Exception: continue
    c=d['chi2']
    if chi is None: chi=numpy.full_like(c, numpy.nan)
    m=numpy.isfinite(c)&~numpy.isfinite(chi); chi[m]=c[m]
print("%d %d" % (int(numpy.isfinite(chi).sum()), chi.size))
PYEOF
}

echo "============================================================"
echo "  STUBBORN finisher (approach C, LEDA)"
echo "  folder : $(pwd)"
echo "  script : $PY"
echo "  up to $MAX_PASSES scan passes until grid is complete"
echo "============================================================"

pass=0
while :; do
    read -r FIN TOT <<< "$(grid_status)"
    echo ""
    echo ">>> grid status before pass $((pass+1)): ${FIN}/${TOT} cells finite"
    if [ "$TOT" -gt 0 ] && [ "$FIN" -eq "$TOT" ]; then
        echo ">>> grid is COMPLETE -- no more scanning needed."
        break
    fi
    if [ "$pass" -ge "$MAX_PASSES" ]; then
        echo "!!! reached MAX_PASSES=$MAX_PASSES with ${FIN}/${TOT} done."
        echo "!!! grid still incomplete -- NOT finalizing. Inspect logs/chunkC*_LEDA.log"
        exit 2
    fi
    pass=$((pass+1))
    echo ">>> scan pass $pass : launching $NCHUNKS workers (auto-resume skips done cells)"

    export OMP_NUM_THREADS=$THREADS_PER_CHUNK
    declare -a PIDS
    for i in $(seq 0 $((NCHUNKS-1))); do
        LOG="logs/chunkC${i}_LEDA.log"
        python -u "$PY" phase=grid chunk=$i nchunks=$NCHUNKS >> "$LOG" 2>&1 &
        PIDS[$i]=$!
    done
    # wait for all; we don't abort on per-chunk failure -- the loop condition
    # (grid completeness) decides whether to go again.
    for i in $(seq 0 $((NCHUNKS-1))); do
        wait "${PIDS[$i]}" || echo "    chunk $i exited non-zero (will be re-checked by grid status)"
    done

    echo ">>> pass $pass done; settling ${SLEEP_BETWEEN}s before re-check"
    sleep "$SLEEP_BETWEEN"
done

# ============================================================
#  postgrid -- merge + verify(@chi2min) + final(@prior).
#  The python script itself asserts the grid is complete and refuses
#  to finalize on a hole, so this is double-guarded.
# ============================================================
echo ""
echo ">>> postgrid: merge + verify + final (FINAL_BY=${FINAL_BY:-prior})"
export OMP_NUM_THREADS=$THREADS_TOTAL
python -u "$PY" phase=postgrid nchunks=$NCHUNKS ${FINAL_BY:+FINAL_BY=$FINAL_BY} 2>&1 \
    | tee logs/postgridC_LEDA.log

echo ""
echo "============================================================"
echo "  LEDA approach-C FINISHED in $(pwd)"
echo "  grid_chi2_GH_C1.npz  -- complete chi^2(v_halo) curve"
echo "  final model saved at the cosmological prior (FINAL_BY=prior)"
echo "============================================================"
