"""
Орбитные диагностики из траекторий: анизотропия beta(r) и intrinsic edge-on
проекции, окрашенные по mass-weighted возрасту/металличности.

  1. {GALAXY}_anisotropy_profiles.pdf — beta_sph(r), beta_z(r), beta_phi(r),
     overall + по орбитальным компонентам (Cappellari+2007 Fig.14, Zhu+2018b
     Fig.8). beta_z = 1 - sigma_z^2/sigma_R^2 ; beta_phi = 1 - sigma_phi^2/sigma_R^2 ;
     beta_sph = 1 - (sigma_theta^2 + sigma_phi^2)/(2 sigma_r^2).
  2. {GALAXY}_edgeon_age.pdf / _edgeon_met.pdf — intrinsic edge-on (x,z) карты,
     цвет = средний по массе возраст / [M/H] (Poci+2019 Fig.8-9, Ding+2023 Fig.5).

Источник данных — те же объекты, что и в chemo_unified.py: скрипт импортирует
chemo_unified и переиспользует архив орбит, веса лучшей модели, разбивку на
ячейки (R, λz) и кэш model_age/model_met. Число орбит N берётся из длины массива
весов (НЕ хардкодится), путь к орбитам — cu.ORBITS_DIR. Поэтому при переходе на
финальную модель (20k орбит) достаточно обновить конфиг chemo_unified.py
(npz_file, orbits_dir) — скрипт подхватит новое N автоматически.

schwarzlib НЕ используется (читаются только готовые траектории + веса).

ВНИМАНИЕ по стоимости: скрипт читает по одному файлу на орбиту с НЕнулевым весом.
Для валидации используйте --stride k (брать каждую k-ю орбиту). Полный прогон —
минуты (зависит от числа орбит с ненулевым весом).

Запуск из папки галактики:
    python orbit_diagnostics.py --stride 10      # быстрая проверка
    python orbit_diagnostics.py                  # полный прогон
"""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


COMP_COLORS = {"corotating": "tab:blue", "spherical": "tab:green",
               "counterrotating": "tab:red", "all": "k"}
COMP_ORDER = ["all", "corotating", "spherical", "counterrotating"]
# Границы λz симметрично ±0.35 (как ст.1/wirb): counter < -0.35 ; spheroid [-0.35,+0.35] ; corotating > +0.35
LZ_LO, LZ_HI = -0.35, 0.35


def read_orbit(path):
    """Быстрое чтение orbit_*.txt -> (Npoints, 6): x,y,z,vx,vy,vz."""
    with open(path) as f:
        a = np.fromstring(f.read(), sep=" ")
    return a.reshape(-1, 6)


def component_of(lz):
    if lz < LZ_LO:
        return "counterrotating"
    if lz < LZ_HI:
        return "spherical"
    return "corotating"


def per_orbit_chem(cu, model_age, model_met):
    """Возраст/металличность каждой орбиты через её ячейку (R, λz)."""
    Rmean = np.asarray(cu.Rmean_list, float)
    lz = np.asarray(cu.lambda_z_list, float)
    idx_r = np.round((Rmean / cu.R_max) * 20).astype(int)
    idx_lz = np.round((1.0 + lz) * 10).astype(int)   # round, как chemo_unified/schwarzlib
    n = len(lz)
    age = np.full(n, np.nan)
    met = np.full(n, np.nan)
    ok = (idx_r >= 0) & (idx_r <= 20) & (idx_lz >= 0) & (idx_lz <= 20)
    ci = np.full(n, -1, dtype=int)
    ci[ok] = cu.dyn_comps_data_map[idx_r[ok], idx_lz[ok]]
    has = ci >= 0
    age[has] = model_age[ci[has]]
    met[has] = model_met[ci[has]]
    return age, met


def accumulate(cu, orb_idx, weights, lz_all, age_all, met_all,
               r_edges, img_half, img_n):
    """Один проход по орбитам: суммы для beta(r) и edge-on карт."""
    nb = len(r_edges) - 1
    comps = ["all", "corotating", "spherical", "counterrotating"]
    # суммы скоростных моментов по радиусу: на компонент
    S = {c: {k: np.zeros(nb) for k in
             ("w", "vR", "vR2", "vphi", "vphi2", "vz", "vz2",
              "vr", "vr2", "vth", "vth2")} for c in comps}
    # edge-on гистограммы (x,z)
    Hw   = {c: np.zeros((img_n, img_n)) for c in comps}
    Hage = {c: np.zeros((img_n, img_n)) for c in comps}
    Hmet = {c: np.zeros((img_n, img_n)) for c in comps}
    img_edges = np.linspace(-img_half, img_half, img_n + 1)

    n_used = 0
    for oi in orb_idx:
        w = float(weights[oi])
        if w <= 0:
            continue
        path = os.path.join(cu.ORBITS_DIR, f"orbit_{oi}.txt")
        if not os.path.exists(path):
            continue
        try:
            p = read_orbit(path)
        except Exception:
            continue
        if p.shape[0] == 0:
            continue
        x, y, z, vx, vy, vz = p.T
        Rcyl = np.hypot(x, y)
        Rcyl = np.where(Rcyl < 1e-6, 1e-6, Rcyl)
        r = np.sqrt(x * x + y * y + z * z)
        r = np.where(r < 1e-6, 1e-6, r)

        vR = (x * vx + y * vy) / Rcyl
        vphi = (x * vy - y * vx) / Rcyl
        vr = (x * vx + y * vy + z * vz) / r
        vth = (z * (x * vx + y * vy) / Rcyl - Rcyl * vz) / r

        comp = component_of(lz_all[oi])
        bidx = np.digitize(r, r_edges) - 1
        m = (bidx >= 0) & (bidx < nb)
        bidx = bidx[m]
        ww = np.full(bidx.shape, w)
        for c in ("all", comp):
            Sc = S[c]
            Sc["w"]   += np.bincount(bidx, ww,            minlength=nb)
            Sc["vR"]  += np.bincount(bidx, ww * vR[m],    minlength=nb)
            Sc["vR2"] += np.bincount(bidx, ww * vR[m]**2, minlength=nb)
            Sc["vphi"]  += np.bincount(bidx, ww * vphi[m],    minlength=nb)
            Sc["vphi2"] += np.bincount(bidx, ww * vphi[m]**2, minlength=nb)
            Sc["vz"]  += np.bincount(bidx, ww * vz[m],    minlength=nb)
            Sc["vz2"] += np.bincount(bidx, ww * vz[m]**2, minlength=nb)
            Sc["vr"]  += np.bincount(bidx, ww * vr[m],    minlength=nb)
            Sc["vr2"] += np.bincount(bidx, ww * vr[m]**2, minlength=nb)
            Sc["vth"]  += np.bincount(bidx, ww * vth[m],    minlength=nb)
            Sc["vth2"] += np.bincount(bidx, ww * vth[m]**2, minlength=nb)

        # edge-on (x, z)
        hw, _, _ = np.histogram2d(x, z, bins=[img_edges, img_edges], weights=np.full_like(x, w))
        a_oi, m_oi = age_all[oi], met_all[oi]
        for c in ("all", comp):
            Hw[c] += hw
            if np.isfinite(a_oi):
                Hage[c] += hw * a_oi
            if np.isfinite(m_oi):
                Hmet[c] += hw * m_oi
        n_used += 1

    return S, Hw, Hage, Hmet, img_edges, n_used


def sigma(S, key2, key1):
    """sigma из сумм: sqrt(<v^2> - <v>^2), безопасно."""
    w = S["w"]
    good = w > 0
    mean = np.where(good, S[key1] / np.where(good, w, 1), 0.0)
    m2 = np.where(good, S[key2] / np.where(good, w, 1), 0.0)
    var = np.clip(m2 - mean**2, 0, None)
    return np.sqrt(var), good


def plot_anisotropy(cu, S, r_cent, outfile):
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), sharex=True)
    defs = [
        (r"$\beta = 1 - (\sigma_\theta^2+\sigma_\phi^2)/(2\sigma_r^2)$", "sph"),
        (r"$\beta_z = 1 - \sigma_z^2/\sigma_R^2$", "z"),
        (r"$\beta_\phi = 1 - \sigma_\phi^2/\sigma_R^2$", "phi"),
    ]
    for ax, (title, kind) in zip(axes, defs):
        for c in COMP_ORDER:
            Sc = S[c]
            sR, gR = sigma(Sc, "vR2", "vR")
            if kind == "z":
                sN, gN = sigma(Sc, "vz2", "vz")
                beta = 1 - (sN**2) / np.where(sR > 0, sR**2, np.nan)
            elif kind == "phi":
                sN, gN = sigma(Sc, "vphi2", "vphi")
                beta = 1 - (sN**2) / np.where(sR > 0, sR**2, np.nan)
            else:
                sr, _ = sigma(Sc, "vr2", "vr")
                sth, _ = sigma(Sc, "vth2", "vth")
                sph, _ = sigma(Sc, "vphi2", "vphi")
                beta = 1 - (sth**2 + sph**2) / np.where(sr > 0, 2 * sr**2, np.nan)
            good = Sc["w"] > 0
            ax.plot(r_cent[good], beta[good], color=COMP_COLORS[c],
                    lw=(2.0 if c == "all" else 1.3),
                    label=cu.COMP_LABELS_RU.get(c, c) if c != "all" else "Все")
        ax.axhline(0, color="0.6", lw=0.8, ls=":")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r"$r$ (arcsec)")
        ax.set_ylim(-1.5, 1.05)
    axes[0].set_ylabel(r"анизотропия $\beta$")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle(f"{cu.GALAXY}: профили анизотропии скоростей")
    fig.tight_layout()
    fig.savefig(outfile, format="pdf")
    plt.close(fig)
    print(f"  сохранено: {outfile}")


def plot_edgeon(cu, Hw, Hval, img_edges, vmin, vmax, value_label, outfile, cmap):
    ext = [img_edges[0], img_edges[-1], img_edges[0], img_edges[-1]]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 9.0))
    for ax, c in zip(axes.ravel(), COMP_ORDER):
        w = Hw[c]
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(w > 0, Hval[c] / w, np.nan)
        # порог по массе: убрать почти пустые пиксели
        thr = np.nanpercentile(w[w > 0], 40) if np.any(w > 0) else 0
        mean = np.where(w >= thr, mean, np.nan)
        im = ax.imshow(np.ma.masked_invalid(mean.T), origin="lower", extent=ext,
                       cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
        ttl = "Все" if c == "all" else cu.COMP_LABELS_RU.get(c, c)
        ax.set_title(ttl, color=COMP_COLORS[c], fontsize=10)
        ax.set_xlabel("x (arcsec)")
        ax.set_ylabel("z (arcsec)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label(value_label)
    fig.suptitle(f"{cu.GALAXY}: intrinsic edge-on — {value_label}", y=1.0)
    fig.tight_layout()
    fig.savefig(outfile, format="pdf")
    plt.close(fig)
    print(f"  сохранено: {outfile}")


def main():
    ap = argparse.ArgumentParser(description="Орбитные диагностики: beta(r) + edge-on.")
    ap.add_argument("--stride", type=int, default=1,
                    help="Брать каждую k-ю орбиту (для быстрой валидации).")
    ap.add_argument("--nbins", type=int, default=14, help="Число радиальных бинов beta.")
    ap.add_argument("--img", type=int, default=60, help="Размер edge-on карты в пикселях.")
    args = ap.parse_args()

    import chemo_unified as cu

    lsq = cu._rload("lsq_solution")
    if lsq is None:
        raise SystemExit("Нет кэша lsq_solution — сначала запустите chemo_unified.py.")
    model_age = np.asarray(lsq["model_age"], float)
    model_met = np.asarray(lsq["model_met"], float)
    mf = cu.compute_mass_fractions()
    cu.COMP_LABELS_RU.update(cu.make_comp_labels_ru(mf))

    weights = np.asarray(cu.weights[cu.model_index], float)
    lz_all = np.asarray(cu.lambda_z_list, float)
    N = len(weights)
    age_all, met_all = per_orbit_chem(cu, model_age, model_met)

    nonzero = int(np.sum(weights > 0))
    orb_idx = np.arange(0, N, args.stride)
    print(f"N орбит={N}, с ненулевым весом={nonzero}, stride={args.stride} "
          f"-> кандидатов к чтению={len(orb_idx)}")

    # радиальные бины и размер edge-on по характерному радиусу модели
    Rchar = float(cu.R_max)
    r_edges = np.linspace(0, 1.5 * Rchar, args.nbins + 1)
    r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
    img_half = 1.3 * Rchar

    S, Hw, Hage, Hmet, img_edges, n_used = accumulate(
        cu, orb_idx, weights, lz_all, age_all, met_all,
        r_edges, img_half, args.img)
    print(f"Прочитано орбит: {n_used}")

    plot_anisotropy(cu, S, r_cent, f"{cu.GALAXY}_anisotropy_profiles.pdf")
    plot_edgeon(cu, Hw, Hage, img_edges, cu.age_min, cu.age_max,
                "возраст (млрд лет)", f"{cu.GALAXY}_edgeon_age.pdf", cmap="jet")
    plot_edgeon(cu, Hw, Hmet, img_edges, cu.met_min, cu.met_max,
                "[M/H] (dex)", f"{cu.GALAXY}_edgeon_met.pdf", cmap="jet")


if __name__ == "__main__":
    main()
