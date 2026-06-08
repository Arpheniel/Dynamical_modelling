"""
⚠️⚠️⚠️ УСТАРЕЛО — ВЫВОД ЭТОГО СКРИПТА НЕВЕРЕН, НЕ ИСПОЛЬЗОВАТЬ ⚠️⚠️⚠️
Ошибка: берёт chem[6]=MAP_VEL1 за «газ». На деле MAP_VEL1 = бо́льшая часть
ЗВЁЗДНОГО света (corr с DAP STELLAR_VEL=+0.99); настоящий газ Hα (DAP
EMLINE_GVEL) = MAP_VEL2 и противовращается. Правильная проверка —
`verify_gas_vs_stars_manga.py`. Итог: газ с ПРОТИВОВРАЩАЮЩИМСЯ диском. См.
INTERPRETATION_STATE.md §2. Оставлен только как запись о том, что было сделано.

Газовый вопрос: с каким диском вращается Hα-газ.

Газ привязан к LOSVD_ID=1 (проверено: все эмиссии, вкл. Hα) => кинематика
газа = MAP_VEL1 (chem[6]). Вопрос: VEL1 совпадает с модельной КОРОТАЦИЕЙ
(λz>0, массивный диск) или КОНТРВРАЩЕНИЕМ (λz<0)?

Транспонирование .T поворачивает дипольное поле скорости на 90°:
  - неверная ориентация  -> |corr| ~ 0 (поля ортогональны),
  - верная ориентация    -> |corr| велик; знак (+corot / -counter) = ответ.
Это даёт критерий, не зависящий от соглашения об ориентации.

Запускать из рабочей папки галактики (pgc_final_C / leda_final_C), где лежит
chemo_unified.py с нужным GALAXY и относительными путями к данным.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import chemo_unified as cu

GAL = cu.GALAXY
print(f"\n===== ГАЗОВЫЙ ВОПРОС: {GAL} =====")

# ---- наблюдательные карты (frame fits, без .T) ----
VEL1 = np.array(cu.chem[6].data, dtype=float)   # gas/LOSVD1
VEL2 = np.array(cu.chem[14].data, dtype=float)  # вторая звёздная компонента
SIG1 = np.array(cu.chem[8].data, dtype=float)

# нормировка скоростей на своё среднее по валидным пикселям (системная)
def demean(a):
    a = a.copy()
    m = np.isfinite(a)
    a[m] -= np.nanmean(a[m])
    return a

VEL1d = demean(VEL1)
VEL2d = demean(VEL2)

# ---- модельные карты ----
corot   = cu.DYN_COMP_LOSVD_MAP_INDICES(14, 21)   # λz>0
counter = cu.DYN_COMP_LOSVD_MAP_INDICES(0, 7)     # λz<0
allcomp = cu.DYN_COMP_LOSVD_MAP_INDICES(0, 21)

mV_corot   = np.array(corot[0],   dtype=float)
mV_counter = np.array(counter[0], dtype=float)
mV_all     = np.array(allcomp[0], dtype=float)
mS_all     = np.array(allcomp[1], dtype=float)


def corr(a, b):
    """Pearson r по общим конечным пикселям."""
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return np.nan, int(m.sum())
    x = a[m] - a[m].mean()
    y = b[m] - b[m].mean()
    den = np.sqrt((x * x).sum() * (y * y).sum())
    if den == 0:
        return np.nan, int(m.sum())
    return float((x * y).sum() / den), int(m.sum())


print("\nКорреляция наблюдаемой скорости с МОДЕЛЬНОЙ (обе ориентации):")
print(f"{'':16s} {'corot asis':>11s} {'corot .T':>11s} {'counter asis':>13s} {'counter .T':>11s}   N")
for name, obs in [("VEL1 (ГАЗ)", VEL1d), ("VEL2", VEL2d)]:
    rca, n  = corr(obs, mV_corot)
    rct, _  = corr(obs, mV_corot.T)
    rua, _  = corr(obs, mV_counter)
    rut, _  = corr(obs, mV_counter.T)
    print(f"{name:16s} {rca:11.3f} {rct:11.3f} {rua:13.3f} {rut:11.3f}  {n}")

# σ как независимый якорь регистрации (знако-чётная величина)
rsa, _ = corr(SIG1, mS_all)
rst, _ = corr(SIG1, mS_all.T)
print(f"\nσ-якорь:  corr(SIG1, model σ_all)  asis={rsa:.3f}  .T={rst:.3f}"
      f"   (большая |corr| указывает верную регистрацию)")

# ---- визуальная сетка ----
panels = [
    ("obs VEL1 (ГАЗ)",     VEL1d,            "vel"),
    ("obs VEL2",           VEL2d,            "vel"),
    ("mod corot  asis",    mV_corot,         "vel"),
    ("mod corot  .T",      mV_corot.T,       "vel"),
    ("mod counter asis",   mV_counter,       "vel"),
    ("mod counter .T",     mV_counter.T,     "vel"),
    ("mod all   asis",     mV_all,           "vel"),
    ("mod all   .T",       mV_all.T,         "vel"),
]
allv = np.concatenate([np.asarray(p[1]).ravel() for p in panels])
allv = allv[np.isfinite(allv)]
vlim = np.nanpercentile(np.abs(allv), 98)

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
for ax, (title, data, kind) in zip(axes.ravel(), panels):
    im = ax.imshow(data, origin="lower", cmap="RdBu_r", vmin=-vlim, vmax=vlim)
    ax.set_title(title, fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
fig.suptitle(f"{GAL}: газ(VEL1) vs модельные диски — поиск верной ориентации .T",
             fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = f"gas_disk_check_{GAL}.png"
fig.savefig(out, dpi=110)
print(f"\nКарта сохранена: {out}")
