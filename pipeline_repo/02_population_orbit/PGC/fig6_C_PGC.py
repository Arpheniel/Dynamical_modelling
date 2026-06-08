"""
fig6_pgc.py — построение mass_dest (распределение массы по радиусу,
разделённое по λ_z компонентам) и Rad_vs_Lambda_z (2D-распределение
орбит в плоскости R_mean × λ_z) для модели PGC 35706.

ИСПРАВЛЕНИЯ относительно исходной версии (см. также fig6_leda_fixed.py):
  1. Убраны магические множители (*7 при формировании orbits[i][1] и /10
     в накоплении масс). Масса в orbits[i][1] теперь просто weights[i] *
     M_stars_total / 1e8 — то есть в единицах 10^8 M_sun. Единицы оси Y на
     mass_dest исправлены на 10^8 M_sun.
  2. Алгоритм binning'а в orb_weights_dist полностью переписан: один
     цикл по орбитам, для каждой определяется бин (i_R, i_λz), и в
     него прибавляется масса. Не O(N × steps²), а O(N).
  3. Раньше орбита могла попасть в несколько бинов (или ни в один) из-за
     перекрытия границ — теперь каждая орбита идёт ровно в один бин.
  4. ax.get_xticks() для секундарной оси заменено на secondary_xaxis с
     functions= аргументами (matplotlib сам пересчитает тики).
  5. Переименовано: mass_polar_ring → mass_spherical (соответствие label'у).
  6. Унифицировано с LEDA: всегда делю archive["Rmean"] на R_scale (в архиве
     Rmean хранится в кпк — см. schwarzlib.py:964).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker as tcr

# ===== КОНСТАНТЫ И ПАРАМЕТРЫ =====
M_stars_total = 5.10e10  # [M_sun] финал C: L_total(4.98e9) × Υ(10.245)
                         # = Υ × densityStars.totalMass(), вычислено через
                         # compute_dynamical_mass.py (Υ=12.4, MGE+distance из run-script).
                         # ВНИМАНИЕ: это масса всей галактики, не только в IFU.
                         # Сумма орбит в IFU даст M_stars_IFU ≈ 5.3×10^10 M_sun
                         # (см. ниже), что согласуется с табличной M_dyn=5.5×10^10.
Max_rad = 10  # [угл. сек.]

# Размер MaNGA IFU (радиус шестиугольника, в угл. сек.) — определяем
# автоматически из bin_scheme файла. Это физически правильный предел
# показа модели: за ним нет наблюдательных данных, и сравнение модели
# с наблюдениями невозможно.
try:
    _bin_scheme = np.loadtxt("bins_PGC35706_Damirs.txt")
    _r_bin = np.sqrt(_bin_scheme[:, 0]**2 + _bin_scheme[:, 1]**2)
    IFU_RADIUS_ARCSEC = _r_bin.max() * 1.0
    print(f"Радиус MaNGA IFU = {IFU_RADIUS_ARCSEC:.2f} угл. сек. "
          f"(из bins_PGC35706_Damirs.txt)")
except (FileNotFoundError, OSError):
    # Fallback если bin_scheme не доступен
    IFU_RADIUS_ARCSEC = Max_rad
    print(f"[WARN] bins_PGC35706_Damirs.txt не найден, использую Max_rad={Max_rad}")

# Используем IFU_RADIUS_ARCSEC как фактический предел графиков.
# Округляем вверх до целого, чтобы тики на оси выглядели аккуратно.
PLOT_R_MAX = int(np.ceil(IFU_RADIUS_ARCSEC))
Filename = "M1e+07_O0_Rh71.6_Vh162_i42_a0_N20000_R1.00_GH_DensitySphHarm.npz"
model_index = 5
distance = 117490  # [kpc]
R_scale = distance * np.pi / 648000  # 1 угл. сек. → R_scale кпк

# Параметры halo и BH для расчёта M_dyn в IFU (финал C, подход C).
RHALO    = 71.6    # r_s NFW, arcsec (космолог. приор)
VHALO    = 161.9   # V_max halo, км/с (best-fit грида)
HALOTYPE = "NFW"   # подход C — NFW
M_BH     = 1e7     # масса BH, M_sun
UPSILON  = 10.245  # M/L (уже учтён в M_stars_total)

# Порог |λ_z| для разделения дисковых vs сферических орбит.
# Должен совпадать с разбиением, использованным в wirb_clean.py.
LAMBDA_Z_DISK_THRESHOLD = 0.35

# ===== ЗАГРУЗКА ДАННЫХ =====
archive = np.load(Filename, allow_pickle=True, encoding='latin1')
weights  = archive["weights"][model_index]
Rmean    = archive["Rmean"]
Lambda_z = archive["MOD_Lambda_z"]

# archive["Rmean"] хранится в кпк (это видно в schwarzlib.py:964 — Rm считается
# напрямую из t[:,0:3] в Agama-единицах, которые kpc). Чтобы отобразить в
# угл. сек., делим на R_scale.
print(f"archive['Rmean'] (кпк): min={Rmean.min():.3f}, max={Rmean.max():.3f}, "
      f"median={np.median(Rmean):.3f}")
Rmean = Rmean / R_scale
print(f"После /R_scale={R_scale:.3f} (угл. сек): "
      f"min={Rmean.min():.3f}, max={Rmean.max():.3f}, median={np.median(Rmean):.3f}")
print(f"Ожидаемый диапазон: [0, {PLOT_R_MAX}] arcsec\n")

# ===== ВЫЧИСЛЕНИЕ M_stars_IFU (звёздная масса в пределах IFU) =====
# weights нормированы на 1 (sum(weights) = 1, см. schwarzlib.py:1010
# "totalMass = 1.0  # weights are normalized to total mass of unity").
# Поэтому weights[i] * M_stars_total = масса орбиты в M_sun.
#
# M_stars_total в шапке — полная звёздная масса галактики.
# Но в IFU попадают только орбиты с Rmean <= IFU_RADIUS_ARCSEC, и их
# суммарная масса — это масса галактики в пределах MaNGA-апертуры (M_stars_IFU).
# Именно она и будет показываться на графике.
# Точно те же границы, что у последнего бина гистограммы:
# bin_idx = round(R / step_size), in_range = (0 <= bin_idx < steps).
# Это даёт R от -step_size/2 до PLOT_R_MAX + step_size/2 — иначе сумма
# столбиков будет систематически расходиться с M_stars_IFU.
_steps = 21  # должно совпадать с steps в построении гистограммы ниже
_step_size = PLOT_R_MAX / (_steps - 1)
_bin_idx_full = np.round(Rmean / _step_size).astype(int)
_within_ifu = (_bin_idx_full >= 0) & (_bin_idx_full < _steps)
M_stars_IFU = M_stars_total * weights[_within_ifu].sum()
print(f"M_stars_total (вся галактика):           {M_stars_total/1e8:>7.2f} × 10^8 M_sun")
print(f"M_stars_IFU   (звёзды в R<{PLOT_R_MAX} arcsec):  {M_stars_IFU/1e8:>7.2f} × 10^8 M_sun"
      f" (доля {weights[_within_ifu].sum()*100:.1f}% массы)")
print()

# ===== ВЫЧИСЛЕНИЕ M_dyn_IFU (полная динамическая масса в пределах IFU) =====
# Используем agama для построения halo-плотности и считаем enclosed mass.
# Если agama недоступен — fallback: M_DM_IFU = 0 (как первое приближение).
try:
    import agama
    # ВАЖНО: те же units, что в run_forstand_grid_*.py!
    # agama.setUnits(mass=1, length=arcsec2kpc, velocity=1)
    arcsec2kpc_units = distance * np.pi / 648000  # = R_scale
    agama.setUnits(mass=1, length=arcsec2kpc_units, velocity=1)

    if HALOTYPE.upper() == "LOG":
        densityHalo = agama.schwarzlib.makeDensityLogHalo(RHALO, VHALO)
    elif HALOTYPE.upper() == "NFW":
        densityHalo = agama.schwarzlib.makeDensityNFWHalo(RHALO, VHALO)
    else:
        raise ValueError(f"Unknown HALOTYPE: {HALOTYPE}")

    # enclosedMass(R) принимает R в length-units (= arcsec в нашей системе).
    # До правой границы последнего бина (а не до PLOT_R_MAX) — для согласованности
    _r_DM_outer = PLOT_R_MAX + _step_size / 2
    M_DM_IFU = densityHalo.enclosedMass(_r_DM_outer)
    print(f"M_DM_IFU      (DM halo в R<{PLOT_R_MAX} arcsec):         "
          f"{M_DM_IFU/1e8:>7.2f} × 10^8 M_sun")
except (ImportError, AttributeError) as _e:
    print(f"[WARN] agama недоступен ({_e}); M_DM_IFU=0 (приближение)")
    M_DM_IFU = 0.0

M_dyn_IFU = M_stars_IFU + M_DM_IFU + M_BH
print(f"M_BH          (массивная BH):                     "
      f"{M_BH/1e8:>7.4f} × 10^8 M_sun")
print(f"M_dyn_IFU     (звёзды + DM + BH в R<IFU):     "
      f"{M_dyn_IFU/1e8:>7.2f} × 10^8 M_sun")
print()

# Подготовка данных орбит: оставляем только те, чей вес значимый.
# orbits[i] = (R_mean_arcsec, mass_1e8_Msun, lambda_z)
WEIGHT_THRESHOLD = 1e-5
orbits = []
for i in range(len(weights)):
    if weights[i] >= WEIGHT_THRESHOLD:
        orbits.append((Rmean[i], weights[i] * M_stars_total / 1e8, Lambda_z[i]))
orbits = np.array(orbits, dtype=float)
print(f"Загружено {len(orbits)} орбит с весом >= {WEIGHT_THRESHOLD}")
print(f"Сумма масс орбит: {orbits[:, 1].sum():.2f} × 10^8 M_sun "
      f"(M_stars_total = {M_stars_total/1e8:.1f} × 10^8 M_sun)")

# ===== ГРАФИК mass_dest: масса по компонентам vs радиус =====
steps = 21
rad_centres = np.linspace(0, PLOT_R_MAX, steps)
zone = PLOT_R_MAX / steps / 2  # полуширина бина по радиусу

# Векторизованное накопление: для каждой орбиты определяем R-бин один раз.
step_size = PLOT_R_MAX / (steps - 1) if steps > 1 else 1.0
R = orbits[:, 0]
M = orbits[:, 1]
Lz = orbits[:, 2]

bin_idx = np.round(R / step_size).astype(int)
in_range = (bin_idx >= 0) & (bin_idx < steps)
n_outside = (~in_range).sum()
if n_outside > 0:
    print(f"  [info] {n_outside} орбит выпадают за пределы IFU [0, {PLOT_R_MAX}] arcsec.")

mass_co_rot       = np.zeros(steps)  # λ_z >= +threshold
mass_counter_rot  = np.zeros(steps)  # λ_z <= -threshold
mass_spherical    = np.zeros(steps)  # |λ_z| < threshold
mass_total        = np.zeros(steps)  # все орбиты

for k in range(steps):
    in_bin = in_range & (bin_idx == k)
    if not np.any(in_bin):
        continue
    M_bin  = M[in_bin]
    Lz_bin = Lz[in_bin]
    mass_total[k]       = M_bin.sum()
    mass_co_rot[k]      = M_bin[Lz_bin >=  LAMBDA_Z_DISK_THRESHOLD].sum()
    mass_counter_rot[k] = M_bin[Lz_bin <= -LAMBDA_Z_DISK_THRESHOLD].sum()
    mass_spherical[k]   = M_bin[np.abs(Lz_bin) < LAMBDA_Z_DISK_THRESHOLD].sum()


# Построение
fig, ax = plt.subplots()
ax.step(rad_centres, mass_spherical,   where='mid',
        label="сферический компонент",          c="grey",  lw=2, ls="--")
ax.step(rad_centres, mass_co_rot,      where='mid',
        label="совращающийся дисковый компонент", c="blue", lw=2, ls="--")
ax.step(rad_centres, mass_counter_rot, where='mid',
        label="противовращающийся дисковый компонент", c="red", lw=2, ls="--")
ax.step(rad_centres, mass_total,       where='mid',
        label="сумма компонентов",              c="black", lw=2, ls="--")

ax.set_xlabel(r"$R_{mean}$, угл. сек.", fontsize=14)
ax.set_ylabel(r"Масса, $10^{8}\ M_{\odot}$", fontsize=14)
ax.set_xlim(0, PLOT_R_MAX)
# Расширяем ylim, чтобы легенда (4 строки) не перекрывала пики гистограммы.
# Множитель 1.6 даёт ~40% свободного пространства сверху.
ax.set_ylim(0, max(mass_total.max() * 1.6, 1e-3))

ax.annotate('PGC 35706', (0.95, 0.95), xycoords='axes fraction',
            fontsize=10, style="italic", ha='right')
ax.tick_params(direction="in", which="both", labelsize=10)
ax.grid(True, alpha=0.3)

# Верхняя ось — кпк
ax_top = ax.secondary_xaxis('top', functions=(lambda x: x * R_scale,
                                              lambda x: x / R_scale))
ax_top.set_xlabel(r'$R_{mean}$, кпк', fontsize=14)
ax_top.tick_params(labelsize=10)

ax.legend(loc='upper left', fontsize=10, frameon=False)

fig.savefig("mass_dest_PGC35706.pdf", format="pdf", bbox_inches='tight')
plt.close(fig)
print("Сохранено mass_dest_PGC35706.pdf")


# ===== ГРАФИК Rad_vs_Lambda_z: 2D-распределение орбит =====
def orb_weights_dist(orbits, R_max, lambda_z_min, lambda_z_max, steps,
                     plot_filename, galaxy_name):
    """Строит 2D-карту массы орбит в плоскости (R_mean, λ_z).

    orbits: массив N×3 — (R_mean_arcsec, mass_1e8_Msun, lambda_z).
    R_max, lambda_z_min, lambda_z_max — границы карты.
    steps — число ячеек по каждой оси.
    """
    R  = orbits[:, 0]
    M  = orbits[:, 1]
    Lz = orbits[:, 2]

    # Шаг по каждой оси.
    R_step  = R_max / (steps - 1)
    Lz_step = (lambda_z_max - lambda_z_min) / (steps - 1)

    # Индексы бинов.
    i_R  = np.round(R / R_step).astype(int)
    i_Lz = np.round((Lz - lambda_z_min) / Lz_step).astype(int)

    # В пределах карты.
    in_map = (i_R >= 0) & (i_R < steps) & (i_Lz >= 0) & (i_Lz < steps)

    image = np.zeros((steps, steps))
    # Для imshow с origin='upper' (default) первая ось = вертикаль сверху вниз.
    # Хотим: λ_z=+1 сверху, λ_z=-1 снизу → инвертируем индекс по Y.
    image_y_idx = (steps - 1) - i_Lz[in_map]
    image_x_idx = i_R[in_map]
    np.add.at(image, (image_y_idx, image_x_idx), M[in_map])

    # Построение
    fig, ax = plt.subplots(figsize=(8, 6))
    pc = ax.imshow(image, cmap="inferno", aspect='auto')

    # Форматтеры тиков: пиксельные индексы → физические значения.
    def x_idx_to_arcsec(x_idx, pos):
        return f'{x_idx * R_step:.0f}'

    def y_idx_to_lambda_z(y_idx, pos):
        # y_idx инвертирован: y_idx=0 → λ_z=+1, y_idx=steps-1 → λ_z=-1.
        return f'{lambda_z_max - y_idx * Lz_step:.1f}'

    ax.xaxis.set_major_formatter(tcr.FuncFormatter(x_idx_to_arcsec))
    ax.yaxis.set_major_formatter(tcr.FuncFormatter(y_idx_to_lambda_z))
    ax.xaxis.set_minor_formatter(tcr.NullFormatter())
    ax.yaxis.set_major_locator(tcr.IndexLocator(2, offset=0.4))
    ax.xaxis.set_minor_locator(tcr.IndexLocator(1, offset=0.5))

    ax.set_xlabel(r"$R_{mean}$, угл. сек.", fontsize=14)
    ax.set_ylabel(r"Циркулярность, $\lambda_{z}$", fontsize=14)
    ax.annotate(galaxy_name, (0.05, 0.95), xycoords='axes fraction',
                fontsize=12, style="italic", verticalalignment='top')

    # Верхняя ось — кпк (через secondary_xaxis с functions= аргументами:
    # тогда matplotlib сам пересчитает тики).
    ax_top = ax.secondary_xaxis(
        'top',
        functions=(lambda x_idx: x_idx * R_step * R_scale,
                   lambda kpc:   kpc / R_scale / R_step))
    ax_top.set_xlabel(r'$R_{mean}$, кпк', fontsize=14)
    ax_top.tick_params(labelsize=10)

    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label(r'Масса, $10^{8}\ M_{\odot}$', fontsize=14)
    cbar.ax.tick_params(labelsize=10)

    ax.tick_params(labelsize=10)
    fig.savefig(plot_filename, format="pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"Сохранено {plot_filename}")


orb_weights_dist(orbits,
                 R_max=PLOT_R_MAX,
                 lambda_z_min=-1.0, lambda_z_max=+1.0,
                 steps=21,
                 plot_filename="Rad_vs_Lambda_z_PGC35706.pdf",
                 galaxy_name="PGC 35706")
