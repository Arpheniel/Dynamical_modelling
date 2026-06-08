# show_status_C.sh -- progress monitor for an approach-C run (PGC or LEDA).
#
# Run it FROM THE SAME folder as forstand_C_*.py and the *_C1 npz files.
# It auto-detects the galaxy from the script name, so the one file works in
# both the PGC and the LEDA working directory.
#
#   bash show_status_C.sh                 # one snapshot
#   watch -n 30 'bash show_status_C.sh'   # refresh every 30 s (if watch exists)
#   TAIL=6 bash show_status_C.sh          # show 6 last log lines per chunk
#
# What it reports
#   - which chunk workers are alive, and whether each is mid-cell right now
#   - per-chunk finite-cell count and the tail of its scan log (crash markers
#     CRASHED / Traceback / launch-failure are highlighted)
#   - the merged v_halo scan table (chi2, best Upsilon) with best point so far
#   - how long since the last cell actually completed (stale-run detector)

set -u

PY=${PY:-python}
SCRIPT=$(ls forstand_C_*.py 2>/dev/null | head -1)
if [ -z "$SCRIPT" ]; then
    echo "no forstand_C_*.py in $(pwd)"; exit 1
fi
GAL=$(echo "$SCRIPT" | sed -E 's/^forstand_C_([A-Za-z0-9]+)_.*/\1/')
NCHUNKS=${NCHUNKS:-4}
TAIL=${TAIL:-3}

echo "================ approach-C status : $GAL   ($(date '+%Y-%m-%d %H:%M:%S')) ================"

# ---------- workers ----------
echo
echo "-- chunk workers --"
any_alive=0
for c in $(seq 0 $((NCHUNKS-1))); do
    ppid=$(pgrep -f "$SCRIPT phase=grid chunk=$c " 2>/dev/null | head -1)
    cpid=$(pgrep -f "phase=onecell .* chunk=$c " 2>/dev/null | head -1)
    log="logs/chunkC${c}_${GAL}.log"
    if [ -n "$ppid" ]; then
        any_alive=1
        if [ -n "$cpid" ]; then state="ALIVE pid=$ppid (computing a cell, child=$cpid)";
        else                    state="ALIVE pid=$ppid (between cells)"; fi
    else
        state="not running"
    fi
    # crash markers in the log
    if [ -f "$log" ]; then
        nbad=$(grep -c -E 'CRASHED|Traceback|could not launch|after .* crashed attempts' "$log" 2>/dev/null)
        marks=""; [ "${nbad:-0}" -gt 0 ] && marks="  [!! $nbad crash/error marks]"
    else
        marks="  [no log yet]"
    fi
    echo "  chunk $c: $state$marks"
    if [ -f "$log" ]; then
        tail -n "$TAIL" "$log" | sed 's/^/      | /'
    fi
done
[ "$any_alive" -eq 0 ] && echo "  (no live workers -- scan finished, crashed, or not started)"

# ---------- merged grid table ----------
echo
echo "-- v_halo scan (merged over master + all chunk npz) --"
"$PY" - <<'PYEOF'
import glob, time, numpy
files = sorted(set(glob.glob('grid_chi2_*_C1_chunk*.npz') + glob.glob('grid_chi2_*_C1.npz')))
if not files:
    print("  no grid_chi2_*_C1*.npz yet -- nothing computed so far"); raise SystemExit
rhalo = vhalo = None
chi2 = ups = None
last_t = 0.0
for f in files:
    try: d = numpy.load(f)
    except Exception as e:
        print("  [warn] cannot read %s (%s)" % (f, e)); continue
    if rhalo is None:
        rhalo, vhalo = d['rhalo'], d['vhalo']
        chi2 = numpy.full((len(rhalo), len(vhalo)), numpy.nan)
        ups  = numpy.full_like(chi2, numpy.nan)
    if d['rhalo'].shape != rhalo.shape or d['vhalo'].shape != vhalo.shape:
        print("  [warn] %s has a different grid shape; skipped" % f); continue
    ci = d['chi2']
    ui = d['upsilon'] if 'upsilon' in d.files else numpy.full_like(ci, numpy.nan)
    m = numpy.isfinite(ci) & ~numpy.isfinite(chi2)
    chi2[m] = ci[m]; ups[m] = ui[m]
    if 'last_progress_t' in d.files:
        last_t = max(last_t, float(d['last_progress_t']))

rs = float(rhalo[0]) if len(rhalo) == 1 else float('nan')
done = int(numpy.isfinite(chi2).sum()); tot = chi2.size
print("  fixed r_s = %.2f arcsec    points done: %d / %d" % (rs, done, tot))
print("  %-10s %-12s %-8s" % ("v_halo", "chi2", "best Y"))
print("  " + "-"*32)
flat_v = vhalo
order = numpy.argsort(flat_v)
col = 0  # single r_s row -> column over vhalo
for j in order:
    c = chi2[0, j]; u = ups[0, j]
    cs = "%.2f" % c if numpy.isfinite(c) else "   --   "
    usv = "%.2f" % u if numpy.isfinite(u) else "  -- "
    print("  %-10.1f %-12s %-8s" % (flat_v[j], cs, usv))
if done:
    j_b = int(numpy.nanargmin(chi2[0]))
    print("  " + "-"*32)
    print("  BEST: v_halo=%.1f km/s   chi2=%.3f   Y=%.3f" %
          (flat_v[j_b], chi2[0, j_b], ups[0, j_b]))
if last_t > 0:
    age = time.time() - last_t
    h = int(age // 3600); m = int((age % 3600) // 60)
    print("  last completed cell: %dh%02dm ago" % (h, m) +
          ("   <-- looks stalled, consider relaunch" if age > 3600 else ""))
PYEOF

echo
echo "(relaunch a stalled/partial scan with:  bash launch_C_${GAL}.sh  -- done points are skipped)"
