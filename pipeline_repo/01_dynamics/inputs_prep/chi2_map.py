"""
chi2_map.py — карта chi^2 динамического грид-скана по гало-параметрам.

Строит из grid_chi2_GH.npz (rhalo x vhalo) две панели:
  - Delta chi^2 = chi^2 - chi^2_min с уровнем доверия 1 sigma = sqrt(2 N_obs)
    (Jin+2024 Fig.4, Lipka & Thomas 2021 Fig.10);
  - карта Upsilon (M*/L) по тем же ячейкам.

N_obs (число кинематических ограничений) = N_bins x N_moments определяется
из kinem-файла и числа GH-моментов (nmom из имени chem-файла). При желании
переопределяется через --n-obs.

schwarzlib НЕ используется — берётся только готовый grid_chi2_GH.npz, поэтому
скрипт корректно отрабатывает локально. После финального грида достаточно
обновить grid_chi2_GH.npz и перезапустить — имена/оси не меняются.

Запуск из папки галактики (PGC / LEDA):
    python chi2_map.py
    python chi2_map.py --grid grid_chi2_GH.npz --n-obs 2704
"""

import os
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import ticker as tcr


def detect_n_moments():
    """Число GH-моментов из имени chem-файла results_*_nmomN.fits (по умолч. 4)."""
    for f in glob.glob("results_*_nmom*.fits"):
        base = os.path.basename(f)
        i = base.find("nmom")
        if i >= 0:
            digits = ""
            for ch in base[i + 4:]:
                if ch.isdigit():
                    digits += ch
                else:
                    break
            if digits:
                return int(digits)
    return 4


def detect_n_obs(n_moments):
    """N_obs = N_kin_bins x N_moments. N_kin_bins — число строк kinem-файла."""
    cands = glob.glob("kinem_gh_*.txt")
    if not cands:
        return None, None
    # предпочесть файл с 'Damir' (рабочая кинематика), иначе первый
    kin = next((c for c in cands if "Damir" in c), cands[0])
    n_bins = np.loadtxt(kin).shape[0]
    return n_bins * n_moments, kin


def make_chi2_map(grid_file, n_obs, galaxy_tag, outfile):
    d = np.load(grid_file, allow_pickle=True)
    rhalo = np.asarray(d["rhalo"], float)
    vhalo = np.asarray(d["vhalo"], float)
    chi2 = np.asarray(d["chi2"], float)            # (n_rhalo, n_vhalo)
    upsilon = np.asarray(d["upsilon"], float)
    incl = float(d["incl"]) if "incl" in d.files else None

    chi2_min = np.nanmin(chi2)
    i_best, j_best = np.unravel_index(np.nanargmin(chi2), chi2.shape)
    dchi2 = chi2 - chi2_min

    # Уровни доверия (Schwarzschild-стандарт: Delta chi^2 = n * sqrt(2 N_obs))
    sig1 = np.sqrt(2.0 * n_obs) if n_obs else None

    # сетка для pcolormesh (по краям)
    def edges(x):
        x = np.asarray(x, float)
        e = np.empty(len(x) + 1)
        e[1:-1] = 0.5 * (x[:-1] + x[1:])
        e[0] = x[0] - (x[1] - x[0]) / 2
        e[-1] = x[-1] + (x[-1] - x[-2]) / 2
        return e

    Rx, Vy = edges(rhalo), edges(vhalo)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.2))

    # ── Левая панель: Delta chi^2 ────────────────────────────────────────────
    vmax = (8.0 * sig1) if sig1 else np.nanpercentile(dchi2, 60)
    pcm = axL.pcolormesh(Rx, Vy, dchi2.T, cmap="viridis_r",
                         vmin=0, vmax=vmax, shading="flat")
    cb = fig.colorbar(pcm, ax=axL)
    cb.set_label(r"$\Delta\chi^2 = \chi^2 - \chi^2_\mathrm{min}$")

    # контуры уровней доверия
    if sig1:
        levels = [sig1, 2 * sig1, 3 * sig1]
        cs = axL.contour(rhalo, vhalo, dchi2.T, levels=levels,
                         colors="white", linewidths=[1.8, 1.2, 0.8])
        fmt = {levels[0]: r"$1\sigma$", levels[1]: r"$2\sigma$", levels[2]: r"$3\sigma$"}
        axL.clabel(cs, fmt=fmt, fontsize=8)

    axL.plot(rhalo[i_best], vhalo[j_best], marker="*", color="red", ms=16,
             mec="k", mew=0.6, zorder=5)
    axL.set_xscale("log")
    axL.set_xlabel(r"$r_\mathrm{halo}$ (arcsec)")
    axL.set_ylabel(r"$v_\mathrm{halo}$ (км/с)")
    title = rf"{galaxy_tag}: $\chi^2_\mathrm{{min}}={chi2_min:.0f}$"
    if sig1:
        title += rf",  $1\sigma:\ \Delta\chi^2={sig1:.1f}$  ($N_\mathrm{{obs}}={n_obs}$)"
    axL.set_title(title, fontsize=9)

    # ── Правая панель: Upsilon ───────────────────────────────────────────────
    pcm2 = axR.pcolormesh(Rx, Vy, upsilon.T, cmap="plasma", shading="flat")
    cb2 = fig.colorbar(pcm2, ax=axR)
    cb2.set_label(r"$\Upsilon$ (M$_\odot$/L$_\odot$)")
    axR.plot(rhalo[i_best], vhalo[j_best], marker="*", color="cyan", ms=16,
             mec="k", mew=0.6, zorder=5)
    axR.set_xscale("log")
    axR.set_xlabel(r"$r_\mathrm{halo}$ (arcsec)")
    axR.set_ylabel(r"$v_\mathrm{halo}$ (км/с)")
    axR.set_title(rf"$\Upsilon$ best-fit $= {upsilon[i_best, j_best]:.2f}$ "
                  rf"при $r_h={rhalo[i_best]:.0f}$, $v_h={vhalo[j_best]:.0f}$", fontsize=9)
    for ax in (axL, axR):
        ax.xaxis.set_major_formatter(tcr.ScalarFormatter())

    fig.tight_layout()
    fig.savefig(outfile, format="pdf")
    plt.close(fig)
    print(f"  сохранено: {outfile}")
    print(f"  best-fit: r_halo={rhalo[i_best]:.1f} arcsec, v_halo={vhalo[j_best]:.0f} км/с, "
          f"Upsilon={upsilon[i_best, j_best]:.2f}, chi2_min={chi2_min:.1f}")


def main():
    ap = argparse.ArgumentParser(description="Карта chi^2 динамического грида.")
    ap.add_argument("--grid", default="grid_chi2_GH.npz")
    ap.add_argument("--n-obs", type=int, default=None,
                    help="Число наблюдательных ограничений (по умолч. авто).")
    ap.add_argument("--galaxy", default=None, help="Тег для заголовка (по умолч. из CWD).")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    galaxy_tag = args.galaxy or os.path.basename(os.getcwd())
    outfile = args.out or f"{galaxy_tag}_chi2_map.pdf"

    nmom = detect_n_moments()
    n_obs = args.n_obs
    if n_obs is None:
        n_obs, kin = detect_n_obs(nmom)
        if n_obs:
            print(f"  N_obs авто: {kin} ({n_obs // nmom} бинов) x {nmom} GH-момента = {n_obs}")
        else:
            print("  ВНИМАНИЕ: kinem-файл не найден, N_obs неизвестно — "
                  "уровни доверия не строятся. Задайте --n-obs.")

    make_chi2_map(args.grid, n_obs, galaxy_tag, outfile)


if __name__ == "__main__":
    main()
