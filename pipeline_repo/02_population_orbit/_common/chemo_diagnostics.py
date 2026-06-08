"""
Стандартные популяционно-орбитальные диагностики для статьи-2.

Скрипт строит набор графиков, которых пока нет среди готовых PDF, но которые
в литературе считаются стандартом для population-orbit decomposition
(Poci+2019, Zhu+2020, Jin+2024 и пр.):

  1. {GALAXY}_circularity_lambda_z.pdf — распределение циркулярности λz
     (дифференциальное + кумулятивное f(<λz)), с массовыми долями компонент.
  2. {GALAXY}_age_metallicity.pdf      — диаграмма возраст–[M/H] по компонентам,
     цвет = масса ячейки.
  3. {GALAXY}_mass_age_met_hist.pdf    — mass-weighted гистограммы возраста и
     [M/H], стопкой по компонентам.
  4. {GALAXY}_age_vs_lambda_z.pdf      — возраст vs λz, цвет = [M/H], размер = масса.

Данные берутся из того же источника, что и chemo_unified.py: скрипт импортирует
chemo_unified (под __main__-guard, поэтому main() не запускается) и переиспользует
уже загруженный архив орбит, разбивку на ячейки (R, λz) и кэш решения BVLS
(results_cache_*/lsq_solution.npz). Орбиты и BVLS НЕ пересчитываются.

Запускать из папки галактики (PGC / LEDA), где лежит соответствующий
chemo_unified.py. После финального пересчёта динамики/декомпозиции достаточно
обновить архив и кэш и снова запустить этот скрипт — имена/оси не меняются.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# Цвета и порядок компонент — совпадают с chemo_unified.py
COMP_COLORS = {"corotating": "tab:blue", "spherical": "tab:green", "counterrotating": "tab:red"}
COMP_ORDER  = ["corotating", "spherical", "counterrotating"]

# Границы λz симметрично ±0.35 (как ст.1/wirb; round-биннинг → idx7/14 = ±0.35):
#   контрвращение λz < -0.35 ; сферический [-0.35,+0.35] ; коротация λz > +0.35
REGION_EDGES = (-0.35, 0.35)
REGIONS = [(-1.0, -0.35, "counterrotating"),
           (-0.35,  0.35, "spherical"),
           ( 0.35,  1.0, "corotating")]


def per_cell_mass(cu):
    """Масса каждой (R, λz)-ячейки = сумма динамических весов орбит в ней.
    Выровнено с порядком cu.dyn_comps_data (как в compute_mass_fractions)."""
    w = np.asarray(cu.weights[cu.model_index], dtype=float)
    return np.array([
        float(np.sum(w[cu.sorted_orbs[ir][ilz]])) if cu.sorted_orbs[ir][ilz] else 0.0
        for ir, ilz in cu.dyn_comps_data
    ])


def cell_component_tags(cu):
    """Тег компонента для каждой ячейки по индексу λz (LZ_REGIONS)."""
    idx_lz = cu.dyn_comps_data[:, 1]
    tags = np.empty(len(idx_lz), dtype=object)
    for lo, hi, tag, _ in cu.LZ_REGIONS:
        tags[(idx_lz >= lo) & (idx_lz < hi)] = tag
    return tags


def cell_lambda_z(cu):
    """Центр λz-бина для каждой ячейки (индекс 0..20 -> λz = idx/10 - 1)."""
    return cu.dyn_comps_data[:, 1] / 10.0 - 1.0


def _shade_regions(ax, mf, colors, annotate_y=None):
    """Подсветить три λz-региона и (опц.) подписать массовые доли."""
    for lo, hi, tag in REGIONS:
        ax.axvspan(lo, hi, color=colors[tag], alpha=0.08, zorder=0)
    for edge in REGION_EDGES:
        ax.axvline(edge, color="0.35", ls="--", lw=0.8, zorder=1)
    if annotate_y is not None:
        for lo, hi, tag in REGIONS:
            ax.text((lo + hi) / 2.0, annotate_y, f"{mf[tag] * 100:.1f}%",
                    ha="center", va="top", color=colors[tag],
                    fontsize=9, fontweight="bold")


# ── 1. Распределение циркулярности λz ───────────────────────────────────────
def plot_circularity(cu, mf, outfile):
    w  = np.asarray(cu.weights[cu.model_index], dtype=float)
    lz = np.asarray(cu.lambda_z_list, dtype=float)
    ok = np.isfinite(lz) & np.isfinite(w) & (w > 0)
    lz, w = lz[ok], w[ok]
    lz = np.clip(lz, -1.0, 1.0)
    wn = w / w.sum()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 7.2), sharex=True)

    bins = np.linspace(-1.0, 1.0, 41)
    ax1.hist(lz, bins=bins, weights=wn, color="0.6", edgecolor="k", linewidth=0.4)
    ax1.set_ylabel(r"Массовая доля на бин $\lambda_z$")
    ax1.set_title(f"{cu.GALAXY}: распределение циркулярности орбит")
    ymax = ax1.get_ylim()[1]
    _shade_regions(ax1, mf, COMP_COLORS, annotate_y=ymax * 0.97)
    ax1.set_ylim(0, ymax * 1.05)

    order = np.argsort(lz)
    cdf = np.cumsum(wn[order])
    ax2.plot(lz[order], cdf, color="k", lw=1.6)
    ax2.set_xlabel(r"$\lambda_z$ (циркулярность орбиты)")
    ax2.set_ylabel(r"$f(<\lambda_z)$")
    ax2.set_xlim(-1.0, 1.0)
    ax2.set_ylim(0, 1.0)
    _shade_regions(ax2, mf, COMP_COLORS)

    fig.tight_layout()
    fig.savefig(outfile, format="pdf")
    plt.close(fig)
    print(f"  сохранено: {outfile}")


# ── 2. Диаграмма возраст–металличность по компонентам ───────────────────────
def plot_age_metallicity(cu, age, met, mass, tags, outfile):
    mpos = mass > 0
    vmin = max(mass[mpos].min(), mass[mpos].max() * 1e-3)
    norm = mcolors.LogNorm(vmin=vmin, vmax=mass[mpos].max())

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3), sharex=True, sharey=True)
    sc = None
    for ax, tag in zip(axes, COMP_ORDER):
        sel = (tags == tag) & mpos
        sc = ax.scatter(met[sel], age[sel], c=mass[sel], norm=norm, cmap="viridis",
                        s=30, edgecolor="k", linewidth=0.25)
        ax.set_title(cu.COMP_LABELS_RU[tag], color=COMP_COLORS[tag], fontsize=10)
        ax.set_xlabel("[M/H] (dex)")
        ax.set_xlim(cu.met_min, cu.met_max)
        ax.set_ylim(cu.age_min, cu.age_max)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Возраст (млрд лет)")
    cbar = fig.colorbar(sc, ax=axes, fraction=0.03, pad=0.02)
    cbar.set_label("Масса ячейки (отн. ед.)")
    fig.suptitle(f"{cu.GALAXY}: возраст–металличность по орбитальным компонентам", y=1.02)
    fig.savefig(outfile, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  сохранено: {outfile}")


# ── 3. Mass-weighted гистограммы возраста и [M/H] ────────────────────────────
def plot_mass_histograms(cu, age, met, mass, tags, outfile):
    total = mass.sum()
    age_bins = np.linspace(cu.age_min, cu.age_max, 16)
    met_bins = np.linspace(cu.met_min, cu.met_max, 16)

    data_age = [age[tags == tag] for tag in COMP_ORDER]
    data_met = [met[tags == tag] for tag in COMP_ORDER]
    w_age = [mass[tags == tag] / total for tag in COMP_ORDER]
    w_met = w_age
    cols = [COMP_COLORS[t] for t in COMP_ORDER]
    labs = [cu.COMP_LABELS_RU[t] for t in COMP_ORDER]

    fig, (axa, axm) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    axa.hist(data_age, bins=age_bins, weights=w_age, stacked=True,
             color=cols, label=labs, edgecolor="k", linewidth=0.3)
    axa.set_xlabel("Возраст (млрд лет)")
    axa.set_ylabel("Массовая доля на бин")
    axa.legend(fontsize=8, loc="best")

    axm.hist(data_met, bins=met_bins, weights=w_met, stacked=True,
             color=cols, edgecolor="k", linewidth=0.3)
    axm.set_xlabel("[M/H] (dex)")
    axm.set_ylabel("Массовая доля на бин")

    fig.suptitle(f"{cu.GALAXY}: распределение массы по возрасту и металличности")
    fig.tight_layout()
    fig.savefig(outfile, format="pdf")
    plt.close(fig)
    print(f"  сохранено: {outfile}")


# ── 4. Возраст vs λz, цвет по [M/H] ──────────────────────────────────────────
def plot_age_vs_lambda(cu, age, met, mass, lz, mf, outfile):
    mpos = mass > 0
    s = 14 + 260.0 * (mass / mass.max())

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    sc = ax.scatter(lz[mpos], age[mpos], c=met[mpos], cmap="coolwarm",
                    s=s[mpos], edgecolor="k", linewidth=0.25,
                    vmin=cu.met_min, vmax=cu.met_max)
    _shade_regions(ax, mf, COMP_COLORS, annotate_y=cu.age_max * 0.98)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(cu.age_min, cu.age_max)
    ax.set_xlabel(r"$\lambda_z$ (циркулярность орбиты)")
    ax.set_ylabel("Возраст (млрд лет)")
    ax.set_title(f"{cu.GALAXY}: возраст vs циркулярность (размер ∝ масса)")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("[M/H] (dex)")
    fig.tight_layout()
    fig.savefig(outfile, format="pdf")
    plt.close(fig)
    print(f"  сохранено: {outfile}")


def main():
    import chemo_unified as cu

    lsq = cu._rload("lsq_solution")
    mf = cu.compute_mass_fractions()
    cu.COMP_LABELS_RU.update(cu.make_comp_labels_ru(mf))

    print("Массовые доли компонент:")
    for tag in COMP_ORDER:
        print(f"  {cu.COMP_LABELS_RU[tag]}")

    # График 1 нужен только архив (веса+λz по орбитам) — строится всегда.
    plot_circularity(cu, mf, f"{cu.GALAXY}_circularity_lambda_z.pdf")

    if lsq is None:
        print("ВНИМАНИЕ: кэш lsq_solution не найден — графики 2–4 (по ячейкам) "
              "пропущены. Сначала запустите chemo_unified.py, чтобы получить "
              "model_age/model_met.")
        return

    age  = np.asarray(lsq["model_age"], dtype=float)
    met  = np.asarray(lsq["model_met"], dtype=float)
    mass = per_cell_mass(cu)
    tags = cell_component_tags(cu)
    lz   = cell_lambda_z(cu)

    assert len(age) == len(mass) == len(tags) == len(lz), \
        f"рассинхрон длин: age={len(age)} mass={len(mass)} tags={len(tags)} lz={len(lz)}"

    plot_age_metallicity(cu, age, met, mass, tags, f"{cu.GALAXY}_age_metallicity.pdf")
    plot_mass_histograms(cu, age, met, mass, tags, f"{cu.GALAXY}_mass_age_met_hist.pdf")
    plot_age_vs_lambda(cu, age, met, mass, lz, mf, f"{cu.GALAXY}_age_vs_lambda_z.pdf")


if __name__ == "__main__":
    main()
