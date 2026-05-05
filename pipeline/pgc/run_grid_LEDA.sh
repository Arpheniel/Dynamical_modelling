#!/bin/bash
# Parallel staged grid launcher for LEDA 2220522.
#
# BEFORE RUNNING: activate the agama virtualenv!
#     agamaenv
#
# LEDA has TWO nearly-degenerate local minima in the adamet chain.
# Strategy: characterise BOTH of them symmetrically with 5x5 grids at
# N=10000, then run N=40000 at each local best-fit so that the orbital
# decomposition maps (paper Fig. 6) can be produced for both solutions
# and compared directly.
#
# Stages:
#   B1  -> 5x5 grid around "upper" minimum (131, 324), step 10" x 15 km/s
#   B2  -> 5x5 grid around "right" minimum (190, 215), step 15" x 20 km/s
#   C1  -> N=40000 with save_orbits at the local best of B1
#   C2  -> N=40000 with save_orbits at the local best of B2
#   ALL -> run B1 + B2 together, then pause for analysis, then run C1/C2
#
# Usage:
#   agamaenv
#   ./run_grid_LEDA.sh ALL              # PARALLEL defaults to 6
#   PARALLEL=12 ./run_grid_LEDA.sh ALL
#   ./run_grid_LEDA.sh C1 135 324       # after B1 identifies best at (135,324)
#   ./run_grid_LEDA.sh C2 190 215
#
# Monitoring from another ssh session:
#   tail -f grid_log_leda.txt
#   ls -1 *_i27_a0_N10000_R0.00_GH_DensityCylindricalLinear.npz | wc -l

set -u
PY=python3
SCRIPT=run_forstand_grid_LEDA_2220522.py
PARALLEL=${PARALLEL:-6}
LOGDIR=logs
MASTER_LOG=grid_log_leda.txt
mkdir -p "$LOGDIR"

if ! python3 -c "import agama" >/dev/null 2>&1; then
    echo "ERROR: 'import agama' failed.  Did you run 'agamaenv' first?" >&2
    exit 1
fi

throttle () {
    while (( $(jobs -rp | wc -l) >= PARALLEL )); do
        wait -n 2>/dev/null || sleep 2
    done
}

submit () {
    local RH=$1 VH=$2 N=$3 SAVE=$4
    local STEM="M1e+07_O0_Rh${RH}_Vh${VH}_i27_a0_N${N}_R0.00_GH_DensityCylindricalLinear"
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

# ----- grid specifications ------
# "upper" minimum: 5x5 around (131, 324), step 10" x 15 km/s
B1_RH=(111 121 131 141 151)
B1_VH=(294 309 324 339 354)

# "right" minimum: 5x5 around (190, 215), step 15" x 20 km/s
# (wider step because the minimum is broader according to the chain)
B2_RH=(160 175 190 205 220)
B2_VH=(175 195 215 235 255)

STAGE=${1:-help}

case "$STAGE" in

help|-h|--help)
    head -30 "$0"
    exit 0
    ;;

B1)
    echo "====== LEDA Stage B1: 5x5 grid around UPPER minimum ($(date)) ======" \
        | tee -a "$MASTER_LOG"
    for RH in "${B1_RH[@]}"; do
        for VH in "${B1_VH[@]}"; do
            submit $RH $VH 10000 no
        done
    done
    wait
    echo "====== Stage B1 done at $(date) ======" | tee -a "$MASTER_LOG"
    ;;

B2)
    echo "====== LEDA Stage B2: 5x5 grid around RIGHT minimum ($(date)) ======" \
        | tee -a "$MASTER_LOG"
    for RH in "${B2_RH[@]}"; do
        for VH in "${B2_VH[@]}"; do
            submit $RH $VH 10000 no
        done
    done
    wait
    echo "====== Stage B2 done at $(date) ======" | tee -a "$MASTER_LOG"
    ;;

ALL)
    echo "====== LEDA ALL: B1 and B2 together ($(date)) ======" \
        | tee -a "$MASTER_LOG"
    for RH in "${B1_RH[@]}"; do
        for VH in "${B1_VH[@]}"; do
            submit $RH $VH 10000 no
        done
    done
    for RH in "${B2_RH[@]}"; do
        for VH in "${B2_VH[@]}"; do
            submit $RH $VH 10000 no
        done
    done
    wait
    echo "====== ALL B1+B2 done at $(date) ======" | tee -a "$MASTER_LOG"
    echo
    echo "NEXT STEP: identify the local minimum in each region from resultsGH.txt"
    echo "and launch C1 / C2 with the chosen (RH, VH):"
    echo "    ./run_grid_LEDA.sh C1 <rh> <vh>     # upper-region best"
    echo "    ./run_grid_LEDA.sh C2 <rh> <vh>     # right-region best"
    ;;

C1)
    FINAL_RH=${2:-131}
    FINAL_VH=${3:-324}
    echo "====== LEDA Stage C1: N=40000 save_orbits at UPPER ($FINAL_RH, $FINAL_VH) ======" \
        | tee -a "$MASTER_LOG"
    submit $FINAL_RH $FINAL_VH 40000 yes
    wait
    ;;

C2)
    FINAL_RH=${2:-190}
    FINAL_VH=${3:-215}
    echo "====== LEDA Stage C2: N=40000 save_orbits at RIGHT ($FINAL_RH, $FINAL_VH) ======" \
        | tee -a "$MASTER_LOG"
    submit $FINAL_RH $FINAL_VH 40000 yes
    wait
    ;;

*)
    echo "Unknown stage: $STAGE.  Run '$0 help' for usage."
    exit 1
    ;;
esac

echo "resultsGH.txt lines: $(wc -l < resultsGH.txt 2>/dev/null || echo 0)" \
    | tee -a "$MASTER_LOG"
