import numpy as np
import matplotlib.pyplot as plt
from matplotlib import transforms, ticker as tcr, colors, patheffects as pe
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
import agama
import scipy
import os
import sys
from tqdm import tqdm
from time import sleep


# ╔══════════════════════════════════════════════════════════════╗
# ║              ВЫБОР ГАЛАКТИКИ (менять здесь)                  ║
# ╚══════════════════════════════════════════════════════════════╝

GALAXY = "PGC"   # "PGC"  — PGC 35706
               # "LEDA" — LEDA 2220522


# ========== ПАРАМЕТРЫ ГАЛАКТИК ==========
# Все значения, отличающиеся между галактиками, собраны здесь.

GALAXY_PARAMS = {
    "PGC": {
        # --- численные параметры модели ---
        "scale":        42,           # размер карты в пикселях
        # ВНИМАНИЕ: "incl" здесь — НЕ наклон галактики. Истинный наклон модели
        # i = 42° (см. имя архива ..._i42_..., forstand_C INCL=42). Это значение —
        # угол поворота для rotation_matrix_x в orbit_map, где конвенция
        # a=0° → edge-on, a=90° → face-on, поэтому a = 90° − i = 90 − 42 = 48°.
        "incl":         48,           # угол поворота orbit_map = 90° − i(=42°); НЕ наклон
        "model_index":  5,            # финал C: Υ=10.245 (best-fit грида), индекс 5 из 7
        "gamma":        -270,         # угол поворота плоскости проекции

        # --- границы BVLS = диапазон наблюд. карт ± малый запас (НЕ края E-MILES) ---
        # данные PGC: age [2.01, 14.00], [M/H] [-0.50, 0.39]
        "age_min":      1.5,
        "age_max":      14.1,
        "met_min":      -0.7,
        "met_max":      0.4,
        "lambda_reg":   0.01,         # финал C: выбрано сканом по χ²_red (колено: railing 82%->24% при +1.5 χ²)

        # --- файлы данных (финальная модель подхода C, 20k орбит) ---
        "npz_file":     "M1e+07_O0_Rh71.6_Vh162_i42_a0_N20000_R1.00_GH_DensitySphHarm.npz",
        "bins_file":    "bins_PGC35706_Damirs.txt",
        "chem_file":    "results_8992-3704_vorb020_md19_ad-1_nmom4.fits",
        "orbits_dir":   "orbits_M1e+07_O0_Rh71.6_Vh162_i42_a0_N2000",  # обрезано до N2000, внутри 20000 орбит

        # --- параметры построения LOSVD-карты ---
        "losvd_shape":      (676, 47),    # (N_bins_spatial, N_vel)
        "bin_scheme_scale": 2,            # множитель координат бин-схемы (x*k + scale/2)
    },
    "LEDA": {
        # --- численные параметры модели ---
        "scale":        32,
        "incl":         63,
        "model_index":  8,
        "gamma":        225,

        # --- границы физических величин ---
        "age_min":      0,
        "age_max":      7,
        "met_min":      -0.75,
        "met_max":      0.5,
        "lambda_reg":   0.05,

        # --- файлы данных ---
        "npz_file":     "M1e+07_O0_Rh131_Vh324_i27_a0_N40000_R0.00_GH_DensityCylindricalLinear.npz",
        "bins_file":    "bins_LEDA_2220522_Damir's.txt",
        "chem_file":    "results_8254-1902_vorb020_md19_ad-1_nmom4.fits",
        "orbits_dir":   "orbits_M1e+07_O0_Rh131_Vh324_i27_a0_N40000",

        # --- параметры построения LOSVD-карты ---
        "losvd_shape":      (263, 47),
        "bin_scheme_scale": 1,            # без умножения: x*1 + scale/2
    },
}

# ---------- Распаковка параметров выбранной галактики ----------
if GALAXY not in GALAXY_PARAMS:
    raise ValueError(f"Неизвестная галактика: '{GALAXY}'. Доступны: {list(GALAXY_PARAMS.keys())}")

P = GALAXY_PARAMS[GALAXY]

scale           = P["scale"]
incl            = P["incl"]   # NB: угол поворота orbit_map (= 90° − i_галактики), НЕ наклон
model_index     = P["model_index"]
gamma           = P["gamma"]
age_min         = P["age_min"]
age_max         = P["age_max"]
met_min         = P["met_min"]
met_max         = P["met_max"]
DEFAULT_LAMBDA_REG  = P["lambda_reg"]
NPZ_FILE        = P["npz_file"]
BINS_FILE       = P["bins_file"]
CHEM_FILE       = P["chem_file"]
ORBITS_DIR      = P["orbits_dir"]
LOSVD_SHAPE     = P["losvd_shape"]
BIN_SCHEME_SCALE = P["bin_scheme_scale"]

print(f"Галактика: {GALAXY}")
print(f"  scale={scale}, incl={incl}, model_index={model_index}, gamma={gamma}")

# ========== КЭШ РЕЗУЛЬТАТОВ ВЫЧИСЛЕНИЙ ==========
# Отдельный кэш для тяжёлых numpy-объектов (weight_matrix, модельные карты, MC).
# Ключ папки включает имя галактики и индекс модели, чтобы не смешивать результаты.
RESULTS_CACHE_DIR = f"results_cache_{GALAXY}_m{model_index}"
if not os.path.exists(RESULTS_CACHE_DIR):
    os.makedirs(RESULTS_CACHE_DIR)


def _rpath(name):
    """Полный путь к файлу кэша результатов."""
    return os.path.join(RESULTS_CACHE_DIR, name + ".npz")


def _rsave(name, **arrays):
    """Сохранить набор массивов в кэш."""
    np.savez_compressed(_rpath(name), **arrays)
    print(f"  [кэш] сохранено: {_rpath(name)}")


def _rload(name):
    """Загрузить кэш; возвращает NpzFile или None если кэша нет."""
    p = _rpath(name)
    if os.path.exists(p):
        print(f"  [кэш] загружено: {p}")
        return np.load(p, allow_pickle=False)
    return None

# ========== ПРОЧИЕ КОНСТАНТЫ И НАСТРОЙКИ ==========
# ── Русские названия компонент (заполняются после вычисления массовых долей) ──
# Структура: { 'tag' | 'All': русское_название }
# Заполняется в main() через make_comp_labels_ru(). До вызова main()
# содержит имена без долей — используется только в функциях, вызываемых из main().
COMP_LABELS_RU: dict = {}


def compute_mass_fractions():
    """
    Вычисляет массовую долю каждого кинематического региона (LZ_REGIONS)
    и полную «All» как сумму.
    Масса компоненты = sum(weights[model_index][orbi]) по всем орбитам ячейки.
    Возвращает dict: {'All': 1.0, tag: fraction, ...}
    """
    idx_lz_all = dyn_comps_data[:, 1]

    # Физическая масса каждого компонента
    comp_mass = np.array([
        np.sum(weights[model_index][sorted_orbs[ir][ilz]])
        if sorted_orbs[ir][ilz] else 0.0
        for ir, ilz in dyn_comps_data
    ])
    total_mass = np.sum(comp_mass)

    fractions = {"All": 1.0}
    for lz_start, lz_end, tag, _ in LZ_REGIONS:
        mask = (idx_lz_all >= lz_start) & (idx_lz_all < lz_end)
        fractions[tag] = np.sum(comp_mass[mask]) / total_mass if total_mass > 0 else 0.0
    return fractions


def make_comp_labels_ru(mass_fractions):
    """
    Строит словарь русских меток для компонент с массовыми долями.
    Вызывать из main() после compute_mass_fractions().
    """
    _base = {
        "All":             "Все компоненты",
        "counterrotating": "Противовращающийся",
        "spherical":       "Сферический",
        "corotating":      "Совращающийся",
    }
    labels = {}
    for key, base_name in _base.items():
        frac = mass_fractions.get(key, None)
        if frac is not None and key != "All":
            labels[key] = f"{base_name} ({frac*100:.1f}%)"
        else:
            labels[key] = base_name
    # Для совместимости с .capitalize() вызовами добавляем и капитализированные ключи
    for tag in list(labels.keys()):
        labels[tag.capitalize()] = labels[tag]
    return labels


# N (число орбит) задаётся ниже, после загрузки архива (= weights.shape[1]).

# Разбивка 21 бина λ_z (индексы 0..20, λ_z = (idx/10) - 1) на три равных региона.
# Каждый регион — 7 бинов; границы выбраны так, чтобы делить [0,21) без остатка.
#   retro   : idx  0.. 6  ->  λ_z ≈ [-1.0, -0.35)   — ретроградные орбиты
#   hot     : idx  7..13  ->  λ_z ≈ [-0.35, +0.35)   — горячие / дисперсионные
#   prograde: idx 14..20  ->  λ_z ≈ [+0.35, +1.0]   — прогрейд/вращение
LZ_REGIONS = [
    (0,  7,  "counterrotating", r"$\lambda_z \in [-1.0,\,-0.35)$"),
    (7,  14, "spherical",       r"$\lambda_z \in [-0.35,\,+0.35)$"),
    (14, 21, "corotating",      r"$\lambda_z \in [+0.35,\,+1.0]$"),
]

# Колормап для химии
cmap_chem = plt.cm.jet
cmap_chem.set_bad(color='white')

# Путь для кэширования карт плотности
CACHE_DIR_ORIGINAL = "density_cache_original_21r"
CACHE_DIR_ROTATED  = "density_cache_rotated_21r"
for cache_dir in [CACHE_DIR_ORIGINAL, CACHE_DIR_ROTATED]:
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)


# ========== ЗАГРУЗКА ДАННЫХ ==========
print("Загрузка данных...")
archive = np.load(NPZ_FILE, allow_pickle=True, encoding='latin1')

weights         = archive["weights"]
N               = weights.shape[1]   # число орбит из архива (финал C: 20000)
ML              = archive["Upsilon"][model_index]
lambda_z_list   = archive["MOD_Lambda_z"]
Rmean_list      = archive["Rmean"]
DYN_COMP_LOSVD  = archive["DYN_COMP_LOSVD"][model_index]

bin_scheme = np.loadtxt(BINS_FILE)
chem       = fits.open(CHEM_FILE)


# ========== ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ==========

def lambda_z_to_index(lambda_z):
    """Преобразование λ_z в индекс (21 бин)"""
    return int(round((1 + lambda_z) * 10))   # round: бины центрированы на целых (как schwarzlib/wirb) → границы idx7/14 = ±0.35


def DYN_COMP_LOSVD_MAP_INDICES(lz_start_idx, lz_end_idx):
    """
    Создание карт кинематики для заданного диапазона индексов λ_z.
    lz_start_idx: начальный индекс (включительно)
    lz_end_idx:   конечный индекс (не включая)
    """
    LOSVD_plt = np.zeros(LOSVD_SHAPE)

    for r in range(21):
        for l_z in range(lz_start_idx, lz_end_idx):
            LOSVD_plt += DYN_COMP_LOSVD[l_z][r]

    GH_moments = agama.ghMoments(
        matrix=LOSVD_plt * ML**-0.5,
        gridv=np.linspace(-250, 250, 46) * ML**0.5,
        degree=2,
        ghorder=6
    )[:, (1, 2, 6, 7, 8, 9)]

    LOSVD_MAP = []

    for gh_moment in range(6):
        KINEM_MAP = np.full((scale, scale), None, dtype=float)
        for bin_i in range(len(bin_scheme)):
            x = int(bin_scheme[bin_i][0] * BIN_SCHEME_SCALE + scale / 2)
            y = int(bin_scheme[bin_i][1] * BIN_SCHEME_SCALE + scale / 2)
            if 0 <= x < scale and 0 <= y < scale:
                KINEM_MAP[y][x] = GH_moments[int(bin_scheme[bin_i][2])][gh_moment]
        LOSVD_MAP.append(KINEM_MAP)

    return LOSVD_MAP


# ========== СОРТИРОВКА ОРБИТ ==========
print("Сортировка орбит...")
sorted_orbs = [[[] for _ in range(21)] for _ in range(21)]   # 21x21 (радиус x λ_z)

dyn_comps_data_map = np.full((21, 21), -1)   # карта radius x lambda_z с индексами компонентов
dyn_comps_data = []
R_max = np.mean(Rmean_list) * 3

for orbi in range(N):
    Rmean    = Rmean_list[orbi]
    lambda_z = lambda_z_list[orbi]
    idx_r    = int(round((Rmean / R_max) * 20))
    idx_lz   = int(round((1 + lambda_z) * 10))   # round (как schwarzlib/wirb): бины на целых, граница ±0.35; floor давал −0.3/+0.4 и рассогласование с DYN_COMP_LOSVD

    if 0 <= idx_r <= 20 and 0 <= idx_lz <= 20 and len(sorted_orbs[idx_r][idx_lz]) == 0:
        dyn_comps_data_map[idx_r][idx_lz] = np.max(dyn_comps_data_map) + 1
        dyn_comps_data.append([idx_r, idx_lz])
        sorted_orbs[idx_r][idx_lz].append(orbi)
    elif 0 <= idx_r <= 20 and 0 <= idx_lz <= 20:
        sorted_orbs[idx_r][idx_lz].append(orbi)

dyn_comps_data = np.array(dyn_comps_data)


# ========== МАТРИЦЫ ПОВОРОТА ==========

def rotation_matrix_x(a):
    return np.array([[1, 0,          0         ],
                     [0, np.cos(a), -np.sin(a)],
                     [0, np.sin(a),  np.cos(a)]])

def rotation_matrix_2d(a):
    """Матрица поворота в 2D плоскости"""
    c = np.cos(np.radians(a))
    s = np.sin(np.radians(a))
    return np.array([[c, -s],
                     [s,  c]])


# ========== КЭШ КАРТ ПЛОТНОСТИ ==========

def get_cache_filename(lambda_z_range, radius_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    lz_min, lz_max = lambda_z_range
    r_min,  r_max  = radius_range
    if plane_rotation == 0:
        return os.path.join(cache_dir, f"density_lz_{lz_min}_{lz_max}_r_{r_min}_{r_max}_incl_{incl}.npy")
    else:
        return os.path.join(cache_dir, f"density_lz_{lz_min}_{lz_max}_r_{r_min}_{r_max}_incl_{incl}_plane_rot_{plane_rotation}.npy")

def save_density_map(density_map, lambda_z_range, radius_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    filename = get_cache_filename(lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    np.save(filename, density_map)
    return filename

def load_density_map(lambda_z_range, radius_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    filename = get_cache_filename(lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    if os.path.exists(filename):
        return np.load(filename)
    return None


# ========== КАРТЫ ПЛОТНОСТИ ==========

def orbit_map_original(orb_groups, incl):
    """Карта плотности БЕЗ вращения плоскости проекции"""
    image = np.zeros((scale, scale))

    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[model_index][orbi]
            try:
                with open(f"{ORBITS_DIR}/orbit_{orbi}.txt", "r") as f:
                    for str_point in f:
                        point = [float(i) for i in str_point.split(" ")]
                        if len(point) < 3:
                            continue
                        position = np.array([point[0], point[1], point[2]])
                        position = position @ rotation_matrix_x(np.radians(incl))
                        x_ = int(position[0] + scale / 2)
                        y_ = int(position[2] + scale / 2)
                        if 0 <= x_ < scale and 0 <= y_ < scale:
                            image[y_][x_] += weight / 1000
            except FileNotFoundError:
                continue

    return image


def orbit_map_rotated_plane(orb_groups, incl, plane_rotation):
    """Карта плотности С вращением плоскости проекции"""
    image      = np.zeros((scale, scale))
    rot_matrix = rotation_matrix_2d(plane_rotation)

    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[model_index][orbi]
            try:
                with open(f"{ORBITS_DIR}/orbit_{orbi}.txt", "r") as f:
                    for str_point in f:
                        point = [float(i) for i in str_point.split(" ")]
                        if len(point) < 3:
                            continue
                        position = np.array([point[0], point[1], point[2]])
                        position = position @ rotation_matrix_x(np.radians(incl))

                        proj_coords    = np.array([position[0], position[2]])
                        rotated_coords = rot_matrix @ proj_coords

                        x_ = int(rotated_coords[0] + scale / 2)
                        y_ = int(rotated_coords[1] + scale / 2)
                        if 0 <= x_ < scale and 0 <= y_ < scale:
                            image[y_][x_] += weight / 1000
            except FileNotFoundError:
                continue

    return image


def create_and_cache_density_map(lambda_z_range, radius_range, incl, plane_rotation=0, use_rotated_plane=False):
    """
    Создание или загрузка карты плотности с кэшированием.
    lambda_z_range: диапазон индексов λ_z (например, (13, 21) для со-вращения)
    """
    cache_dir = CACHE_DIR_ROTATED if use_rotated_plane else CACHE_DIR_ORIGINAL

    cached = load_density_map(lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    if cached is not None:
        return cached

    orb_list = [sorted_orbs[radius][lz]
                for lz     in range(*lambda_z_range)
                for radius in range(*radius_range)]

    if use_rotated_plane:
        image = orbit_map_rotated_plane(orb_list, incl, plane_rotation)
    else:
        image = orbit_map_original(orb_list, incl)

    save_density_map(image, lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    return image


# ========== ПОПУЛЯЦИОННО-ДИНАМИЧЕСКОЕ МОДЕЛИРОВАНИЕ ==========

# Подготовка векторов с наблюдаемой химией
bin_map = chem[5].data
age_map = chem[22].data
age_err = chem[23].data
met_map = chem[24].data
met_err = chem[25].data

# Наблюдаемая кинематика — из того же fits-файла (spectr == chem)
velocity_obs = chem[6].data - np.nanmean(chem[6].data)  # вычитаем среднее (системная скорость)
sigma_obs    = chem[8].data

num_bin      = np.max(bin_map) + 1
dyn_comp_num = len(dyn_comps_data)

bin_met = np.zeros(num_bin * 2 - 3)   # вектор наблюдаемой металличности
bin_age = np.zeros(num_bin * 2 - 3)   # вектор наблюдаемых возрастов

# Ошибки по бинам (только пространственная часть, без регуляризации)
bin_age_err = np.zeros(num_bin)
bin_met_err = np.zeros(num_bin)

for x_ in range(scale):
    for y_ in range(scale):
        if bin_map[y_][x_] != -1:
            bin_met[bin_map[y_][x_]] = met_map[y_][x_]
            bin_age[bin_map[y_][x_]] = age_map[y_][x_]
            bin_age_err[bin_map[y_][x_]] = age_err[y_][x_]
            bin_met_err[bin_map[y_][x_]] = met_err[y_][x_]


def reg_matrix_1(dyn_comp_i, bin_i, dyn_comps_data_map, dyn_comps_data, lambda_reg):
    """
    Матрица регуляризации.
    Возвращает  4 если dyn_comp_i == bin_i,
               -1 если компонент соседствует с bin_i на карте radius x lambda_z.
    """
    if dyn_comp_i == bin_i:
        return 4 * lambda_reg
    neighbors = [(+1, 0), (0, +1), (-1, 0), (0, -1)]
    for dr, dl in neighbors:
        try:
            idx = dyn_comps_data_map[dyn_comps_data[bin_i][0] + dr][dyn_comps_data[bin_i][1] + dl]
            if dyn_comp_i == idx:
                return -1 * lambda_reg
        except Exception:
            pass
    return 0


def create_weight_matrix_of_dyn_comps(lambda_reg):
    """
    Создание матрицы N_bins × N_comps (матрица весов) с нормировкой по массе бина
    и добавлением матрицы регуляризации.
    """
    print("=" * 30)
    print("Создание матрицы весов для популяционно-динамического моделирования:")

    weight_matrix_dyn_comps = np.zeros((dyn_comp_num, num_bin * 2 - 3))
    weight_map_binned       = np.zeros((num_bin,))
    dyn_comp_ind = 0

    for t, lambda_z in zip(tqdm(range(21)), range(21)):
        for radius in range(21):
            if len(sorted_orbs[radius][lambda_z]) == 0:
                continue

            dyn_comp_bins_weights = np.zeros((num_bin * 2 - 3,))
            lambda_z_range = (lambda_z, lambda_z + 1)
            radius_range   = (radius,   radius   + 1)

            dyn_comp_weight_map = create_and_cache_density_map(
                lambda_z_range, radius_range, incl,
                plane_rotation=gamma, use_rotated_plane=True
            )

            for x_ in range(scale):
                for y_ in range(scale):
                    if bin_map[y_][x_] != -1:
                        dyn_comp_bins_weights[bin_map[y_][x_]] += dyn_comp_weight_map[y_][x_]
                        weight_map_binned[bin_map[y_][x_]]     += dyn_comp_weight_map[y_][x_]

            weight_matrix_dyn_comps[dyn_comp_ind] = dyn_comp_bins_weights
            dyn_comp_ind += 1

    print("=" * 30)
    print("Сделано.")
    print("Масса в матрице весов от общей массы галактики: ", np.sum(weight_map_binned))
    print("Количество динамических компонентов: ", dyn_comp_num)
    print("Количество бинов: ", num_bin)

    # Нормализация весов
    for bin_i in range(num_bin):
        for dyn_comp_i in range(dyn_comp_num):
            if weight_matrix_dyn_comps[dyn_comp_i][bin_i] != 0:
                weight_matrix_dyn_comps[dyn_comp_i][bin_i] /= np.sum(weight_map_binned[bin_i])

    weight_matrix_dyn_comps = np.transpose(weight_matrix_dyn_comps)

    # Добавление матрицы регуляризации
    for dyn_comp_i in range(dyn_comp_num):
        for bin_i in range(num_bin - 3):
            weight_matrix_dyn_comps[bin_i + num_bin][dyn_comp_i] = reg_matrix_1(
                dyn_comp_i, bin_i, dyn_comps_data_map, dyn_comps_data, lambda_reg
            )

    return weight_matrix_dyn_comps, weight_map_binned


# ========== ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ ГРАФИКОВ ==========

# Перевод пикселя в угловые секунды относительно центра кадра:
#   bin_scheme хранит координаты в arcsec, и они кладутся в пиксели как
#   px = coord_arcsec * BIN_SCHEME_SCALE + scale/2
#   => coord_arcsec = (px - scale/2) / BIN_SCHEME_SCALE
_arcsec_per_pixel = 1.0 / BIN_SCHEME_SCALE


def _arcsec_formatter(px, pos):
    """Форматтер тиков: пиксельный индекс → смещение в ″ от центра (сдвиг -0.5″)."""
    return f"{(px - scale / 2) * _arcsec_per_pixel + 0.5:.1f}"


def apply_spatial_axes(ax):
    """
    Устанавливает на осях подписи Δx (″) / Δy (″) со смещением от центра кадра.
    Применяется ко всем картам в пиксельных координатах (scale × scale).
    """
    fmt = tcr.FuncFormatter(_arcsec_formatter)
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.set_xlabel("угл. сек.")
    ax.set_ylabel("угл. сек.")


def _joint_limits(model_arr, obs_arr, obs_percentile=2):
    """
    Возвращает (vmin, vmax) для совместной цветовой шкалы модели и наблюдений.

    Наблюдения могут содержать выбросы, поэтому их диапазон обрезается
    перцентилями [obs_percentile, 100-obs_percentile] (по умолчанию 2%–98%).
    Модельные значения берутся целиком (nanmin/nanmax) — они уже ограничены
    физическими bounds lsq_linear и выбросов не содержат.

    Итоговые vmin/vmax — объединение обоих диапазонов.
    """
    obs_lo = np.nanpercentile(obs_arr, obs_percentile)
    obs_hi = np.nanpercentile(obs_arr, 100 - obs_percentile)

    model_lo = np.nanmin(model_arr)
    model_hi = np.nanmax(model_arr)

    return min(obs_lo, model_lo), max(obs_hi, model_hi)


def build_model_chem_image(model_chem, weight_matrix):
    """
    Строит пространственную карту химии из вектора компонент и матрицы весов.
    Возвращает массив (scale, scale) с NaN вне бинов.
    """
    image = np.full((scale, scale), np.nan)
    for x_ in range(scale):
        for y_ in range(scale):
            if bin_map[y_][x_] != -1:
                image[y_][x_] = np.sum(weight_matrix[bin_map[y_][x_]] * model_chem)
    return image


def build_model_chem_image_lz_region(model_chem, weight_matrix, lz_start, lz_end):
    """
    Пространственная карта химии только для компонент с idx_lz в [lz_start, lz_end).

    В каждом пикселе вычисляется взвешенное среднее:
        value = sum(W[bin, mask] * chem[mask]) / sum(W[bin, mask])
    где mask отбирает компоненты данного lambda_z-региона.

    Если суммарный вес региона в пикселе равен нулю (регион не представлен),
    пиксель получает NaN.
    """
    mask = (dyn_comps_data[:, 1] >= lz_start) & (dyn_comps_data[:, 1] < lz_end)

    image = np.full((scale, scale), np.nan)
    for x_ in range(scale):
        for y_ in range(scale):
            b = bin_map[y_][x_]
            if b == -1:
                continue
            w_region = weight_matrix[b][mask]
            total_w  = np.sum(w_region)
            if total_w == 0:
                continue
            image[y_][x_] = np.sum(w_region * model_chem[mask]) / total_w
    return image


# ========== ГРАФИКИ ==========

def tickers_radius_formatter(x, pos):
    return f'{(x / 21) * R_max:.1f}'

def tickers_lambda_z_formatter(y, pos):
    return f'{(y - 10) / 10:.1f}'


def make_lambda_z_vs_radius_vs_chem_plot(model_chem, dyn_comps_data, chem_min, chem_max, chem="chem"):
    """
    График λ_z vs radius, раскрашенный по химии.
    Оси: радиус в arcsec (по x), λ_z (по y) — фазовое пространство, не небесная плоскость.
    """
    image = np.full((21, 21), np.nan)
    for dyn_comp, dyn_comp_i in zip(dyn_comps_data, range(len(dyn_comps_data))):
        image[dyn_comp[1]][dyn_comp[0]] = model_chem[dyn_comp_i]

    fig, ax = plt.subplots()
    pc = ax.imshow(image, cmap=cmap_chem, vmin=chem_min, vmax=chem_max)
    fig.colorbar(pc, ax=ax, label=chem)
    ax.set_ylabel(r"$\lambda_z$")
    ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_lambda_z_formatter))
    ax.set_xlabel('Радиус (″)')
    ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_radius_formatter))
    fig.savefig(f"{GALAXY}_lambda_z_vs_radius_vs_{chem}.pdf", format="pdf")


def make_model_chem_plot(image, chem_min, chem_max, chem="chem"):
    """
    Пространственная карта модельной химии.
    image: уже построенный массив (scale, scale) из build_model_chem_image().
    Оси: Δx, Δy в угловых секундах от центра кадра.
    """
    fig, ax = plt.subplots()
    pc = ax.imshow(image, cmap=cmap_chem, vmin=chem_min, vmax=chem_max, origin="lower")
    fig.colorbar(pc, ax=ax, label=chem)
    apply_spatial_axes(ax)
    fig.savefig(f"{GALAXY}_model_{chem}.pdf", format="pdf")


def make_obs_chem_plot(chem_map, chem_min, chem_max, chem="chem"):
    """
    Пространственная карта наблюдаемой химии.
    Оси: Δx, Δy в угловых секундах от центра кадра.
    """
    fig, ax = plt.subplots()
    pc = ax.imshow(chem_map, cmap=cmap_chem, vmin=chem_min, vmax=chem_max, origin="lower")
    fig.colorbar(pc, ax=ax, label=chem)
    apply_spatial_axes(ax)
    fig.savefig(f"{GALAXY}_obs_{chem}.pdf", format="pdf")


def make_model_dencity_plot(weight_map_binned):
    """
    Пространственная карта модельной поверхностной плотности (log-норма).
    Оси: Δx, Δy в угловых секундах от центра кадра.
    """
    image = np.full((scale, scale), np.nan)
    for x_ in range(scale):
        for y_ in range(scale):
            if bin_map[y_][x_] != -1:
                image[y_][x_] = weight_map_binned[bin_map[y_][x_]]

    fig, ax = plt.subplots()
    pc = ax.imshow(image, cmap=cmap_chem, norm="log", origin="lower")
    fig.colorbar(pc, ax=ax, label="Поверхностная плотность")
    apply_spatial_axes(ax)
    fig.savefig(f"{GALAXY}_model_dencity.pdf", format="pdf")


def make_kinem_plot(kinem_map, label, filename, transpose=False):
    """
    Пространственная карта кинематического момента (V, σ, h3, …).
    Цветовая шкала: nanmin/nanmax самой карты.
    Оси: Δx, Δy в угловых секундах от центра кадра, Y растёт снизу вверх.

    transpose=True  — для модельных карт (KINEM_MAP из DYN_COMP_LOSVD_MAP_INDICES):
                      нужно .T из-за специфики LOSVD-проекции.
    transpose=False — для наблюдательных карт (chem[6], chem[8]):
                      fits-данные уже ориентированы правильно.
    """
    vmin = np.nanmin(kinem_map)
    vmax = np.nanmax(kinem_map)
    data = kinem_map.T if transpose else kinem_map

    fig, ax = plt.subplots()
    pc = ax.imshow(data, cmap="RdBu_r", vmin=vmin, vmax=vmax, origin="lower")
    fig.colorbar(pc, ax=ax, label=label)
    apply_spatial_axes(ax)
    fig.savefig(filename, format="pdf")
    plt.close(fig)



def make_summary_grid(
    age_images, met_images,
    vel_maps,  vel_transpose,
    sig_maps,  sig_transpose,
    age_vmin, age_vmax,
    met_vmin, met_vmax,
    filename="summary_grid.pdf",
):
    """
    Сводная таблица карт: 4 колонки x 5 строк, графики вплотную друг к другу.

    Компоновка:
      - Два вложенных GridSpec: верхний (колорбары) и нижний (карты).
      - hspace=0, wspace=0 на обоих — ноль зазора между ячейками.
      - Тики и подписи осей только на внешних краях таблицы.
      - Подписи строк размещаются как ylabel первой колонки.
      - Колорбары прижаты к верхнему краю карт через figure-координаты.
    """
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    import matplotlib.cm     as mplcm
    import matplotlib.colors as mplcolors

    ROW_LABELS = [
        "Наблюдения",
        "Модель",
        COMP_LABELS_RU.get("corotating", "Совращающийся"),
        COMP_LABELS_RU.get("spherical",  "Сферический"),
        COMP_LABELS_RU.get("counterrotating", "Противовращающийся"),
    ]
    COL_LABELS = [
        "Возраст (млрд лет)",
        "Металличность [Fe/H]",
        r"$V$ (км/с)",
        r"$\sigma$ (км/с)",
    ]
    N_ROWS, N_COLS = 5, 4

    # ── Шкалы ─────────────────────────────────────────────────────────────────
    all_vel = np.concatenate([np.asarray(m).ravel() for m in vel_maps])
    all_sig = np.concatenate([np.asarray(m).ravel() for m in sig_maps])

    vel_abs = max(abs(np.nanpercentile(all_vel, 2)), abs(np.nanpercentile(all_vel, 98)))
    vel_vmin, vel_vmax = -vel_abs, vel_abs

    sig_vmin = np.nanpercentile(all_sig, 2)
    sig_vmax = np.nanpercentile(all_sig, 98)

    # ── Три тика с центром в 0 ───────────────────────────────────────────────
    # Пиксельная позиция нуля (где форматтер даёт 0.0):
    #   (px - scale/2) * _arcsec_per_pixel + 0.5 = 0  =>  px = scale/2 - 0.5/BIN_SCHEME_SCALE
    _px_zero = scale / 2 - 0.5 * BIN_SCHEME_SCALE
    # Шаг в arcsec: 60% от полного радиуса, округлённый до ближайшего "красивого" числа
    _half_arcsec = (scale / 2) * _arcsec_per_pixel
    _raw_step    = _half_arcsec * 0.6
    # Округляем до 1, 2, 5, 10, 15, 20 ...
    _nice = [1, 2, 5, 10, 15, 20, 25, 30]
    _step_arcsec = min(_nice, key=lambda v: abs(v - _raw_step))
    _step_px     = _step_arcsec * BIN_SCHEME_SCALE
    _sym3_ticks  = [_px_zero - _step_px, _px_zero, _px_zero + _step_px]

    cmaps    = [cmap_chem,  cmap_chem,  "RdBu_r",  "RdBu_r" ]
    vmins    = [age_vmin,   met_vmin,   vel_vmin,  sig_vmin ]
    vmaxs    = [age_vmax,   met_vmax,   vel_vmax,  sig_vmax ]
    col_data = [age_images, met_images, vel_maps,  sig_maps ]
    col_trans = [
        [False] * N_ROWS,
        [False] * N_ROWS,
        vel_transpose,
        sig_transpose,
    ]

    # ── Фигура и два GridSpec ─────────────────────────────────────────────────
    # Ячейки карт должны быть квадратными (изображения scale×scale).
    # Фиксируем ширину фигуры, вычисляем высоту из условия квадратности ячеек.
    LEFT   = 0.13
    RIGHT  = 0.98
    BOTTOM = 0.05
    CBAR_H_IN = 0.25   # высота полосы колорбаров в дюймах (фиксированная)
    GAP_IN    = 0.0    # зазор между колорбарами и картами в дюймах

    FIG_W = 14.0       # ширина фигуры в дюймах
    map_w_in = FIG_W * (RIGHT - LEFT)             # ширина зоны карт в дюймах
    cell_w   = map_w_in / N_COLS                  # ширина одной ячейки
    map_h_in = cell_w * N_ROWS                    # высота зоны карт (квадратные ячейки)
    top_margin_in = CBAR_H_IN + GAP_IN + 0.1      # место над картами
    bot_margin_in = FIG_W * BOTTOM                # место под картами (пропорционально)
    FIG_H  = map_h_in + top_margin_in + bot_margin_in

    # Пересчитываем в доли фигуры
    BOTTOM_frac  = bot_margin_in / FIG_H
    TOP_MAPS     = BOTTOM_frac + map_h_in / FIG_H
    BOT_CBAR     = TOP_MAPS + GAP_IN / FIG_H
    TOP_CBAR     = BOT_CBAR + CBAR_H_IN / FIG_H

    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # GridSpec для N_ROWS×N_COLS карт — без зазоров, квадратные ячейки
    gs_maps = GridSpec(
        nrows=N_ROWS, ncols=N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOTTOM_frac, top=TOP_MAPS,
        figure=fig,
    )

    # GridSpec для строки колорбаров — прижата к картам снизу
    gs_cbars = GridSpec(
        nrows=1, ncols=N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOT_CBAR, top=TOP_CBAR,
        figure=fig,
    )

    # ── Колорбары ─────────────────────────────────────────────────────────────
    for col in range(N_COLS):
        cbar_ax = fig.add_subplot(gs_cbars[0, col])
        norm = mplcolors.Normalize(vmin=vmins[col], vmax=vmaxs[col])
        sm   = mplcm.ScalarMappable(cmap=cmaps[col], norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        cb.set_label(COL_LABELS[col], labelpad=4, fontsize=18)
        cb.ax.tick_params(labelsize=14)
        cbar_ax.xaxis.set_ticks_position("top")
        cbar_ax.xaxis.set_label_position("top")

    # ── Карты ─────────────────────────────────────────────────────────────────
    axes = {}
    for row in range(N_ROWS):
        for col in range(N_COLS):
            ax = fig.add_subplot(gs_maps[row, col])
            axes[(row, col)] = ax

            data = col_data[col][row]
            if col_trans[col][row]:
                data = data.T

            ax.imshow(
                data,
                cmap=cmaps[col],
                vmin=vmins[col], vmax=vmaxs[col],
                origin="lower",
                aspect="auto",   # ячейки квадратные благодаря figsize
                interpolation="nearest",
            )

            # ── Тики и подписи ────────────────────────────────────────────────
            fmt = tcr.FuncFormatter(_arcsec_formatter)

            # X-ось: подпись и тики только в последней строке
            if row == N_ROWS - 1:
                ax.xaxis.set_major_formatter(fmt)
                ax.xaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="x", labelsize=14, labelbottom=True)
                ax.set_xlabel("угл. сек.", fontsize=16)
            else:
                ax.tick_params(axis="x", labelbottom=False)
                ax.set_xlabel("")

            # Y-ось: тики и форматтер только в первой колонке
            if col == 0:
                ax.yaxis.set_major_formatter(fmt)
                ax.yaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="y", labelsize=14, labelleft=True)
                # Подпись строки + ось Y совмещены
                ax.set_ylabel(ROW_LABELS[row], fontsize=9)   # 9: длинные подписи компонент не лезут друг на друга
            else:
                ax.tick_params(axis="y", labelleft=False)
                ax.set_ylabel("")


    fig.savefig(filename, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Сводная таблица графиков сохранена: {filename}")



def build_density_image_met_subset(weight_matrix, weight_map_binned, model_met, met_lo, met_hi):
    """
    Карта поверхностной плотности только для компонент с металличностью
    в диапазоне [met_lo, met_hi), в абсолютных единицах (не нормированных).

    weight_matrix[b][comp] = нормированная доля компонента в бине b, поэтому
    умножаем на weight_map_binned[b] — абсолютную массу бина. Тогда сумма
    карт по всем подмножествам точно равна полной карте плотности.
    """
    mask = (model_met >= met_lo) & (model_met < met_hi)
    image = np.full((scale, scale), np.nan)
    for x_ in range(scale):
        for y_ in range(scale):
            b = bin_map[y_][x_]
            if b == -1:
                continue
            val = np.sum(weight_matrix[b][mask]) * weight_map_binned[b]
            if val > 0:
                image[y_][x_] = val
    return image


def DYN_COMP_LOSVD_MAP_MET_SUBSET(model_met, met_lo, met_hi):
    """
    Карты кинематики (V, σ, h3…) только для компонент с металличностью
    в диапазоне [met_lo, met_hi).

    Суммируем LOSVD только по тем (idx_lz, idx_r)-ячейкам, чей компонент
    попадает в подмножество по металличности.
    """
    LOSVD_plt = np.zeros(LOSVD_SHAPE)

    for comp_i, (idx_r, idx_lz) in enumerate(dyn_comps_data):
        if met_lo <= model_met[comp_i] < met_hi:
            LOSVD_plt += DYN_COMP_LOSVD[idx_lz][idx_r]

    GH_moments = agama.ghMoments(
        matrix=LOSVD_plt * ML**-0.5,
        gridv=np.linspace(-250, 250, 46) * ML**0.5,
        degree=2,
        ghorder=6,
    )[:, (1, 2, 6, 7, 8, 9)]

    LOSVD_MAP = []
    for gh_moment in range(6):
        KINEM_MAP = np.full((scale, scale), None, dtype=float)
        for bin_i in range(len(bin_scheme)):
            x = int(bin_scheme[bin_i][0] * BIN_SCHEME_SCALE + scale / 2)
            y = int(bin_scheme[bin_i][1] * BIN_SCHEME_SCALE + scale / 2)
            if 0 <= x < scale and 0 <= y < scale:
                KINEM_MAP[y][x] = GH_moments[int(bin_scheme[bin_i][2])][gh_moment]
        LOSVD_MAP.append(KINEM_MAP)

    return LOSVD_MAP


def make_met_grid(
    weight_matrix,
    weight_map_binned,
    model_met,
    met_vmin, met_vmax,
    n_subsets=4,
    filename="met_grid.pdf",
):
    """
    Сетка карт: 3 колонки (плотность | V | σ) × n_subsets строк.

    Строки — равномерные подмножества по металличности от met_vmin до met_vmax.
    Над каждой колонкой — горизонтальный колорбар.
    Компоновка: hspace=0, wspace=0, тики только на внешних краях.
    """
    from matplotlib.gridspec import GridSpec
    import matplotlib.cm     as mplcm
    import matplotlib.colors as mplcolors

    N_ROWS = n_subsets
    N_COLS = 3

    # ── Границы подмножеств по металличности ─────────────────────────────────
    edges = np.linspace(met_vmin, met_vmax, N_ROWS + 1)
    row_labels = [
        f"[Fe/H] ∈ [{edges[i]:.2f}, {edges[i+1]:.2f})"
        for i in range(N_ROWS)
    ]

    # ── Строим данные для каждой строки ──────────────────────────────────────
    dens_maps = []
    vel_maps  = []
    sig_maps  = []

    for i in range(N_ROWS):
        lo, hi = edges[i], edges[i + 1]
        # последний бин — включаем правую границу
        hi_use = hi + 1e-9 if i == N_ROWS - 1 else hi

        dens_maps.append(build_density_image_met_subset(weight_matrix, weight_map_binned, model_met, lo, hi_use))
        kin = DYN_COMP_LOSVD_MAP_MET_SUBSET(model_met, lo, hi_use)
        vel_maps.append(kin[0])
        sig_maps.append(kin[1])

    # ── Общие шкалы ──────────────────────────────────────────────────────────
    all_dens = np.concatenate([m.ravel() for m in dens_maps])
    dens_vmin = np.nanpercentile(all_dens, 2)
    dens_vmax = np.nanpercentile(all_dens, 98)

    all_vel  = np.concatenate([np.asarray(m).ravel() for m in vel_maps])
    vel_abs  = max(abs(np.nanpercentile(all_vel, 2)), abs(np.nanpercentile(all_vel, 98)))
    vel_vmin, vel_vmax = -vel_abs, vel_abs

    all_sig  = np.concatenate([np.asarray(m).ravel() for m in sig_maps])
    sig_vmin = np.nanpercentile(all_sig, 2)
    sig_vmax = np.nanpercentile(all_sig, 98)

    COL_LABELS = ["Поверхн. плотность", r"$V$ (км/с)", r"$\sigma$ (км/с)"]
    cmaps      = [cmap_chem,   "RdBu_r",        "RdBu_r"           ]
    vmins      = [dens_vmin,   vel_vmin,         sig_vmin           ]
    vmaxs      = [dens_vmax,   vel_vmax,         sig_vmax           ]
    col_data   = [dens_maps,   vel_maps,          sig_maps           ]
    # плотность — без .T; кинематика — с .T (модельные карты)
    col_trans  = [[False]*N_ROWS, [True]*N_ROWS, [True]*N_ROWS]

    # ── Компоновка ────────────────────────────────────────────────────────────
    LEFT     = 0.18
    RIGHT    = 0.98
    BOTTOM   = 0.05
    CBAR_H   = 0.035
    GAP      = 0.005
    TOP_MAPS = 0.96 - CBAR_H - GAP
    BOT_CBAR = TOP_MAPS + GAP
    TOP_CBAR = BOT_CBAR + CBAR_H

    fig = plt.figure(figsize=(10, 3 * N_ROWS + 1))

    gs_maps = GridSpec(
        nrows=N_ROWS, ncols=N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOTTOM, top=TOP_MAPS,
        figure=fig,
    )
    gs_cbars = GridSpec(
        nrows=1, ncols=N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOT_CBAR, top=TOP_CBAR,
        figure=fig,
    )

    # ── Колорбары ─────────────────────────────────────────────────────────────
    for col in range(N_COLS):
        cbar_ax = fig.add_subplot(gs_cbars[0, col])
        norm = mplcolors.Normalize(vmin=vmins[col], vmax=vmaxs[col])
        sm   = mplcm.ScalarMappable(cmap=cmaps[col], norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        cb.set_label(COL_LABELS[col], labelpad=4, fontsize=18)
        cb.ax.tick_params(labelsize=14)
        cbar_ax.xaxis.set_ticks_position("top")
        cbar_ax.xaxis.set_label_position("top")

    # ── Карты ─────────────────────────────────────────────────────────────────
    fmt = tcr.FuncFormatter(_arcsec_formatter)

    _px_zero    = scale / 2 - 0.5 * BIN_SCHEME_SCALE
    _half_arcsec = (scale / 2) * _arcsec_per_pixel
    _raw_step    = _half_arcsec * 0.6
    _nice        = [1, 2, 5, 10, 15, 20, 25, 30]
    _step_arcsec = min(_nice, key=lambda v: abs(v - _raw_step))
    _step_px     = _step_arcsec * BIN_SCHEME_SCALE
    _sym3_ticks  = [_px_zero - _step_px, _px_zero, _px_zero + _step_px]

    for row in range(N_ROWS):
        for col in range(N_COLS):
            ax   = fig.add_subplot(gs_maps[row, col])
            data = col_data[col][row]
            if col_trans[col][row]:
                data = data.T
            ax.imshow(
                data,
                cmap=cmaps[col],
                vmin=vmins[col], vmax=vmaxs[col],
                origin="lower",
                aspect="equal",
                interpolation="nearest",
            )

            # X-ось: только нижняя строка
            if row == N_ROWS - 1:
                ax.xaxis.set_major_formatter(fmt)
                ax.xaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="x", labelsize=14, labelbottom=True)
                ax.set_xlabel("угл. сек.", fontsize=16)
            else:
                ax.tick_params(axis="x", labelbottom=False)
                ax.set_xlabel("")

            # Y-ось: только первая колонка
            if col == 0:
                ax.yaxis.set_major_formatter(fmt)
                ax.yaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="y", labelsize=14, labelleft=True)
                ax.set_ylabel(row_labels[row], fontsize=9)
            else:
                ax.tick_params(axis="y", labelleft=False)
                ax.set_ylabel("")

    fig.savefig(filename, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Сетка по металличности сохранена: {filename}")



def build_density_image_age_subset(weight_matrix, weight_map_binned, model_age, age_lo, age_hi):
    """
    Карта поверхностной плотности только для компонент с возрастом
    в диапазоне [age_lo, age_hi), в абсолютных единицах.

    Аналогично build_density_image_met_subset: умножаем на weight_map_binned[b],
    чтобы получить абсолютную массу, а не нормированную долю.
    """
    mask = (model_age >= age_lo) & (model_age < age_hi)
    image = np.full((scale, scale), np.nan)
    for x_ in range(scale):
        for y_ in range(scale):
            b = bin_map[y_][x_]
            if b == -1:
                continue
            val = np.sum(weight_matrix[b][mask]) * weight_map_binned[b]
            if val > 0:
                image[y_][x_] = val
    return image


def DYN_COMP_LOSVD_MAP_AGE_SUBSET(model_age, age_lo, age_hi):
    """
    Карты кинематики (V, σ, …) только для компонент с возрастом
    в диапазоне [age_lo, age_hi).
    """
    LOSVD_plt = np.zeros(LOSVD_SHAPE)

    for comp_i, (idx_r, idx_lz) in enumerate(dyn_comps_data):
        if age_lo <= model_age[comp_i] < age_hi:
            LOSVD_plt += DYN_COMP_LOSVD[idx_lz][idx_r]

    GH_moments = agama.ghMoments(
        matrix=LOSVD_plt * ML**-0.5,
        gridv=np.linspace(-250, 250, 46) * ML**0.5,
        degree=2,
        ghorder=6,
    )[:, (1, 2, 6, 7, 8, 9)]

    LOSVD_MAP = []
    for gh_moment in range(6):
        KINEM_MAP = np.full((scale, scale), None, dtype=float)
        for bin_i in range(len(bin_scheme)):
            x = int(bin_scheme[bin_i][0] * BIN_SCHEME_SCALE + scale / 2)
            y = int(bin_scheme[bin_i][1] * BIN_SCHEME_SCALE + scale / 2)
            if 0 <= x < scale and 0 <= y < scale:
                KINEM_MAP[y][x] = GH_moments[int(bin_scheme[bin_i][2])][gh_moment]
        LOSVD_MAP.append(KINEM_MAP)

    return LOSVD_MAP


def make_age_grid(
    weight_matrix,
    weight_map_binned,
    model_age,
    age_vmin, age_vmax,
    n_subsets=4,
    filename="age_grid.pdf",
):
    """
    Сетка карт: 3 колонки (плотность | V | σ) × n_subsets строк.

    Строки — равномерные подмножества по возрасту от age_vmin до age_vmax.
    Над каждой колонкой — горизонтальный колорбар.
    Компоновка: hspace=0, wspace=0, тики только на внешних краях.
    """
    from matplotlib.gridspec import GridSpec
    import matplotlib.cm     as mplcm
    import matplotlib.colors as mplcolors

    N_ROWS = n_subsets
    N_COLS = 3

    edges = np.linspace(age_vmin, age_vmax, N_ROWS + 1)
    row_labels = [
        f"Возраст ∈ [{edges[i]:.1f}, {edges[i+1]:.1f}) млрд лет"
        for i in range(N_ROWS)
    ]

    dens_maps = []
    vel_maps  = []
    sig_maps  = []

    for i in range(N_ROWS):
        lo, hi = edges[i], edges[i + 1]
        hi_use = hi + 1e-9 if i == N_ROWS - 1 else hi

        dens_maps.append(build_density_image_age_subset(weight_matrix, weight_map_binned, model_age, lo, hi_use))
        kin = DYN_COMP_LOSVD_MAP_AGE_SUBSET(model_age, lo, hi_use)
        vel_maps.append(kin[0])
        sig_maps.append(kin[1])

    all_dens = np.concatenate([m.ravel() for m in dens_maps])
    dens_vmin = np.nanpercentile(all_dens, 2)
    dens_vmax = np.nanpercentile(all_dens, 98)

    all_vel  = np.concatenate([np.asarray(m).ravel() for m in vel_maps])
    vel_abs  = max(abs(np.nanpercentile(all_vel, 2)), abs(np.nanpercentile(all_vel, 98)))
    vel_vmin, vel_vmax = -vel_abs, vel_abs

    all_sig  = np.concatenate([np.asarray(m).ravel() for m in sig_maps])
    sig_vmin = np.nanpercentile(all_sig, 2)
    sig_vmax = np.nanpercentile(all_sig, 98)

    COL_LABELS = ["Поверхн. плотность", r"$V$ (км/с)", r"$\sigma$ (км/с)"]
    cmaps      = [cmap_chem,   "RdBu_r",        "RdBu_r"           ]
    vmins      = [dens_vmin,   vel_vmin,         sig_vmin           ]
    vmaxs      = [dens_vmax,   vel_vmax,         sig_vmax           ]
    col_data   = [dens_maps,   vel_maps,          sig_maps           ]
    col_trans  = [[False]*N_ROWS, [True]*N_ROWS, [True]*N_ROWS]

    LEFT     = 0.18
    RIGHT    = 0.98
    BOTTOM   = 0.05
    CBAR_H   = 0.035
    GAP      = 0.005
    TOP_MAPS = 0.96 - CBAR_H - GAP
    BOT_CBAR = TOP_MAPS + GAP
    TOP_CBAR = BOT_CBAR + CBAR_H

    fig = plt.figure(figsize=(10, 3 * N_ROWS + 1))

    gs_maps = GridSpec(
        nrows=N_ROWS, ncols=N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOTTOM, top=TOP_MAPS,
        figure=fig,
    )
    gs_cbars = GridSpec(
        nrows=1, ncols=N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOT_CBAR, top=TOP_CBAR,
        figure=fig,
    )

    for col in range(N_COLS):
        cbar_ax = fig.add_subplot(gs_cbars[0, col])
        norm = mplcolors.Normalize(vmin=vmins[col], vmax=vmaxs[col])
        sm   = mplcm.ScalarMappable(cmap=cmaps[col], norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        cb.set_label(COL_LABELS[col], labelpad=4, fontsize=18)
        cb.ax.tick_params(labelsize=14)
        cbar_ax.xaxis.set_ticks_position("top")
        cbar_ax.xaxis.set_label_position("top")

    fmt = tcr.FuncFormatter(_arcsec_formatter)
    _px_zero     = scale / 2 - 0.5 * BIN_SCHEME_SCALE
    _half_arcsec = (scale / 2) * _arcsec_per_pixel
    _raw_step    = _half_arcsec * 0.6
    _nice        = [1, 2, 5, 10, 15, 20, 25, 30]
    _step_arcsec = min(_nice, key=lambda v: abs(v - _raw_step))
    _step_px     = _step_arcsec * BIN_SCHEME_SCALE
    _sym3_ticks  = [_px_zero - _step_px, _px_zero, _px_zero + _step_px]

    for row in range(N_ROWS):
        for col in range(N_COLS):
            ax   = fig.add_subplot(gs_maps[row, col])
            data = col_data[col][row]
            if col_trans[col][row]:
                data = data.T
            ax.imshow(
                data,
                cmap=cmaps[col],
                vmin=vmins[col], vmax=vmaxs[col],
                origin="lower",
                aspect="equal",
                interpolation="nearest",
            )

            if row == N_ROWS - 1:
                ax.xaxis.set_major_formatter(fmt)
                ax.xaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="x", labelsize=14, labelbottom=True)
                ax.set_xlabel("угл. сек.", fontsize=16)
            else:
                ax.tick_params(axis="x", labelbottom=False)
                ax.set_xlabel("")

            if col == 0:
                ax.yaxis.set_major_formatter(fmt)
                ax.yaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="y", labelsize=14, labelleft=True)
                ax.set_ylabel(row_labels[row], fontsize=9)
            else:
                ax.tick_params(axis="y", labelleft=False)
                ax.set_ylabel("")

    fig.savefig(filename, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Сетка по возрасту сохранена: {filename}")

def make_radial_chem_profiles(model_age, model_met, weight_matrix, r_max=None, filename="radial_chem_profiles.pdf"):
    """
    Графики радиус–возраст и радиус–металличность для полной модели
    и трёх кинематических компонент из LZ_REGIONS.

    Радиус усреднён по углу: для каждого idx_r берём все компоненты
    данного региона с этим idx_r и считаем взвешенное среднее,
    где вес компоненты = суммарный вклад в пространственные бины.

    Радиус в угловых секундах: r = idx_r / 20 * R_max.
    """
    # Суммарный вес каждой компоненты (только по пространственным бинам, не по регуляризации)
    # Физические веса: сумма orbit weights по орбитам каждой ячейки
    comp_weights = np.array([
        np.sum(weights[model_index][sorted_orbs[ir][ilz]])
        if sorted_orbs[ir][ilz] else 0.0
        for ir, ilz in dyn_comps_data
    ])

    idx_r_all = dyn_comps_data[:, 0]   # 0..20
    idx_lz_all = dyn_comps_data[:, 1]  # 0..20
    _r_scale      = r_max if r_max is not None else R_max
    r_arcsec_bins = np.arange(21) / 20.0 * _r_scale

    # Наборы: (метка, маска по λ_z, цвет, стиль линии)
    subsets = [(COMP_LABELS_RU.get("All", "Все компоненты"), np.ones(len(dyn_comps_data), dtype=bool), "black", "-",  2.5)]
    for lz_start, lz_end, tag, lz_label in LZ_REGIONS:
        mask = (idx_lz_all >= lz_start) & (idx_lz_all < lz_end)
        colors_map = {"corotating": "tab:blue", "spherical": "tab:green", "counterrotating": "tab:red"}
        subsets.append((COMP_LABELS_RU.get(tag, tag.capitalize()), mask, colors_map[tag], "--", 1.8))

    _r_scale = r_max if r_max is not None else R_max
    r_arcsec_bins = np.arange(21) / 20.0 * _r_scale

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    ax_age, ax_met = axes

    for label, lz_mask, color, ls, lw in subsets:
        age_profile = np.full(21, np.nan)
        met_profile = np.full(21, np.nan)

        for ir in range(21):
            r_mask = idx_r_all == ir
            sel = r_mask & lz_mask
            if not np.any(sel):
                continue
            w = comp_weights[sel]
            w_sum = np.sum(w)
            if w_sum == 0:
                continue
            age_profile[ir] = np.sum(w * model_age[sel]) / w_sum
            met_profile[ir] = np.sum(w * model_met[sel]) / w_sum

        valid = ~np.isnan(age_profile)
        ax_age.step(r_arcsec_bins[valid], age_profile[valid],
                    where="mid", label=label, color=color, ls=ls, lw=lw)
        valid = ~np.isnan(met_profile)
        ax_met.step(r_arcsec_bins[valid], met_profile[valid],
                    where="mid", label=label, color=color, ls=ls, lw=lw)

    for ax, ylabel, title in [
        (ax_age, "Возраст (млрд лет)",         "Радиальный профиль возраста"),
        (ax_met, "Металличность [Fe/H]",  "Радиальный профиль металличности"),
    ]:
        ax.set_xlabel(r"$r$ (″)", fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title, fontsize=14)
        ax.tick_params(labelsize=12)
        ax.legend(fontsize=11)
        ax.set_xlim(left=0)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(filename, format="pdf")
    plt.close(fig)
    print(f"Радиальные профили химии сохранены: {filename}")


def make_radial_chem_profiles_stacked(model_age, model_met, weight_matrix,
                                      r_max=None, filename="radial_chem_profiles_stacked.pdf"):
    """
    Радиальные профили возраста и металличности, разнесённые по вертикали.

    Каждая строка — один кинематический компонент (All / Corotating / Spherical /
    Counterrotating). В каждой строке два subplot: возраст (слева) и металличность
    (справа). Оси X (радиус) одинаковы для всех строк и подписаны только снизу.
    Ступенчатый стиль (where="mid").
    """
    comp_weights = np.array([
        np.sum(weights[model_index][sorted_orbs[ir][ilz]])
        if sorted_orbs[ir][ilz] else 0.0
        for ir, ilz in dyn_comps_data
    ])

    idx_r_all  = dyn_comps_data[:, 0]
    idx_lz_all = dyn_comps_data[:, 1]
    _r_scale      = r_max if r_max is not None else R_max
    r_arcsec_bins = np.arange(21) / 20.0 * _r_scale

    subsets = [(COMP_LABELS_RU.get("All", "Все компоненты"), np.ones(len(dyn_comps_data), dtype=bool), "black",    "-",  2.5)]
    colors_map = {"corotating": "tab:blue", "spherical": "tab:green", "counterrotating": "tab:red"}
    for lz_start, lz_end, tag, _ in LZ_REGIONS:
        mask = (idx_lz_all >= lz_start) & (idx_lz_all < lz_end)
        subsets.append((tag.capitalize(), mask, colors_map[tag], "-", 2.0))

    N = len(subsets)
    fig, axes = plt.subplots(
        nrows=N, ncols=2,
        figsize=(12, 3.5 * N),
        sharex=True,
        sharey="col",       # одинаковый Y в каждой колонке
        gridspec_kw={"hspace": 0.08, "wspace": 0.25},
    )

    for row, (label, lz_mask, color, ls, lw) in enumerate(subsets):
        age_profile = np.full(21, np.nan)
        met_profile = np.full(21, np.nan)

        for ir in range(21):
            sel = (idx_r_all == ir) & lz_mask
            if not np.any(sel):
                continue
            w = comp_weights[sel]
            w_sum = np.sum(w)
            if w_sum == 0:
                continue
            age_profile[ir] = np.sum(w * model_age[sel]) / w_sum
            met_profile[ir] = np.sum(w * model_met[sel]) / w_sum

        ax_age = axes[row, 0]
        ax_met = axes[row, 1]

        valid = ~np.isnan(age_profile)
        ax_age.step(r_arcsec_bins[valid], age_profile[valid],
                    where="mid", color=color, ls=ls, lw=lw)

        valid = ~np.isnan(met_profile)
        ax_met.step(r_arcsec_bins[valid], met_profile[valid],
                    where="mid", color=color, ls=ls, lw=lw)

        ax_age.set_ylabel("Возраст (млрд лет)",        fontsize=12)
        ax_met.set_ylabel("Металличность [Fe/H]", fontsize=12)

        for ax in (ax_age, ax_met):
            ax.set_xlim(left=0)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=11)
            if row == N - 1:
                ax.set_xlabel(r"$r$ (″)", fontsize=13)

        # Подпись компоненты — правая ось строки
        ax_r = ax_met.twinx()
        ax_r.set_yticks([])
        ax_r.set_ylabel(label, fontsize=13, rotation=270, labelpad=20, va="bottom")

    # Подписи колонок сверху
    axes[0, 0].set_title("Возраст (млрд лет)",        fontsize=14)
    axes[0, 1].set_title("Металличность [Fe/H]", fontsize=14)

    fig.savefig(filename, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Стековые радиальные профили сохранены: {filename}")


def make_residuals_grid(
    model_age_image, age_map,
    model_met_image, met_map,
    velocity_map,    velocity_obs,
    sigma_map,       sigma_obs,
    age_vmin, age_vmax,
    met_vmin, met_vmax,
    filename="residuals_grid.pdf",
):
    """
    Сетка невязок (модель − наблюдения).

    Колонки : Возраст | Металличность | V | σ
              (над каждой — свой колорбар, одна физическая величина на колонку)
    Строки  : Наблюдения | Модель | Невязка (модель − наблюдения)

    Шкала наблюдений и модели для каждой колонки одинакова.
    Шкала невязки — симметричная вокруг 0 (RdBu_r), своя на каждую колонку.
    """
    from matplotlib.gridspec import GridSpec
    import matplotlib.cm     as mplcm
    import matplotlib.colors as mplcolors

    # Кинематика модели — нужно .T
    vel_model = velocity_map.T
    sig_model = sigma_map.T

    # Невязки
    res_age = model_age_image - age_map

    # Металличность [M/H] = log10(Z/Zsun) — уже логарифмическая величина.
    # Невязка = разность в dex (= log10 отношения метелличностей), как у возраста.
    # (Старый код брал sign(10^model-10^obs)*log10|10^model-10^obs| — это давало
    #  огромный выброс там, где model≈obs: log10 от крошечной линейной разности.)
    res_met = model_met_image - met_map

    res_vel = vel_model - velocity_obs
    res_sig = sig_model - sigma_obs

    # Шкалы obs/model для кинематики (перцентили 2–98 по обоим вместе)
    def joint_lim(a, b):
        ab = np.concatenate([np.asarray(a).ravel(), np.asarray(b).ravel()])
        return np.nanpercentile(ab, 2), np.nanpercentile(ab, 98)

    vel_vmin, vel_vmax = joint_lim(velocity_obs, vel_model)
    vel_abs = max(abs(vel_vmin), abs(vel_vmax))
    vel_vmin, vel_vmax = -vel_abs, vel_abs          # симметрично для V

    sig_vmin, sig_vmax = joint_lim(sigma_obs, sig_model)

    # Единая шкала на колонку: объединяем obs, model и residual
    # Для V используем симметричную шкалу (RdBu_r), для остальных — joint range.
    def unified_lim(arr_list):
        all_vals = np.concatenate([np.asarray(a).ravel() for a in arr_list])
        return np.nanpercentile(all_vals, 2), np.nanpercentile(all_vals, 98)

    age_u_min, age_u_max = unified_lim([age_map, model_age_image, res_age])
    met_u_min, met_u_max = unified_lim([met_map, model_met_image, res_met])

    vel_all = np.concatenate([velocity_obs.ravel(), vel_model.ravel(), res_vel.ravel()])
    vel_u_abs = max(abs(np.nanpercentile(vel_all, 2)), abs(np.nanpercentile(vel_all, 98)))
    vel_u_min, vel_u_max = -vel_u_abs, vel_u_abs   # симметрично для V

    sig_u_min, sig_u_max = unified_lim([sigma_obs, sig_model, res_sig])

    # cols_data: (obs, model, residual, cmap, vmin, vmax)
    cols_data = [
        (age_map,      model_age_image, res_age, cmap_chem, age_u_min, age_u_max),
        (met_map,      model_met_image, res_met, cmap_chem, met_u_min, met_u_max),
        (velocity_obs, vel_model,       res_vel, "RdBu_r",  vel_u_min, vel_u_max),
        (sigma_obs,    sig_model,       res_sig, "RdBu_r",  sig_u_min, sig_u_max),
    ]

    COL_LABELS = ["Возраст (млрд лет)", "Металличность [Fe/H]", r"$V$ (км/с)", r"$\sigma$ (км/с)"]
    ROW_LABELS = ["Наблюдения", "Модель", "Невязка (мод. − набл.)"]

    N_ROWS, N_COLS = 3, 4

    LEFT     = 0.10
    RIGHT    = 0.98
    BOTTOM   = 0.05
    CBAR_H   = 0.030
    GAP      = 0.004
    TOP_MAPS = 0.94 - CBAR_H - GAP
    BOT_CB   = TOP_MAPS + GAP
    TOP_CB   = BOT_CB   + CBAR_H

    fig = plt.figure(figsize=(16, 10))

    gs_maps = GridSpec(
        N_ROWS, N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOTTOM, top=TOP_MAPS,
        figure=fig,
    )
    gs_cb = GridSpec(
        1, N_COLS,
        hspace=0, wspace=0,
        left=LEFT, right=RIGHT,
        bottom=BOT_CB, top=TOP_CB,
        figure=fig,
    )

    fmt = tcr.FuncFormatter(_arcsec_formatter)
    _px_zero     = scale / 2 - 0.5 * BIN_SCHEME_SCALE
    _half_arcsec = (scale / 2) * _arcsec_per_pixel
    _raw_step    = _half_arcsec * 0.6
    _nice        = [1, 2, 5, 10, 15, 20, 25, 30]
    _step_arcsec = min(_nice, key=lambda v: abs(v - _raw_step))
    _step_px     = _step_arcsec * BIN_SCHEME_SCALE
    _sym3_ticks  = [_px_zero - _step_px, _px_zero, _px_zero + _step_px]

    # ── Один колорбар на колонку ──────────────────────────────────────────────
    for ci, (_, _, _, cmap, vmin, vmax) in enumerate(cols_data):
        cax  = fig.add_subplot(gs_cb[0, ci])
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
        sm   = mplcm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb   = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(COL_LABELS[ci], labelpad=4, fontsize=15)
        cb.ax.tick_params(labelsize=12)
        cax.xaxis.set_ticks_position("top")
        cax.xaxis.set_label_position("top")

    # ── Карты ─────────────────────────────────────────────────────────────────
    for col, (obs, model, res, cmap, vmin, vmax) in enumerate(cols_data):
        row_items = [
            (obs,   cmap, vmin, vmax),
            (model, cmap, vmin, vmax),
            (res,   cmap, vmin, vmax),
        ]
        for row, (data, cmap_use, vmin_use, vmax_use) in enumerate(row_items):
            ax = fig.add_subplot(gs_maps[row, col])
            ax.imshow(
                data,
                cmap=cmap_use, vmin=vmin_use, vmax=vmax_use,
                origin="lower", aspect="equal", interpolation="nearest",
            )

            # X-ось: только нижняя строка
            if row == N_ROWS - 1:
                ax.xaxis.set_major_formatter(fmt)
                ax.xaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="x", labelsize=13, labelbottom=True)
                ax.set_xlabel("угл. сек.", fontsize=14)
            else:
                ax.tick_params(axis="x", labelbottom=False)
                ax.set_xlabel("")

            # Y-ось: только первая колонка
            if col == 0:
                ax.yaxis.set_major_formatter(fmt)
                ax.yaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="y", labelsize=13, labelleft=True)
                ax.set_ylabel(ROW_LABELS[row], fontsize=13)
            else:
                ax.tick_params(axis="y", labelleft=False)
                ax.set_ylabel("")

    fig.savefig(filename, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Карты невязок сохранены: {filename}")


# ── Монте-Карло: оценка устойчивости модели ───────────────────────────────────

def _mc_perturb_age_vec(bin_age_obs, bin_age_err_obs, n_sim):
    """
    Возмущение вектора возрастов (лог-нормальное, E[X] = age_obs, без сдвига Дженсена).

    Формулы:
        sigma_log = sqrt(ln(1 + (err/age)^2))
        mu_log    = ln(age) - 0.5 * sigma_log^2   =>  E[exp(N(mu,s))] = age
    Возвращает массив (n_sim, num_bin).
    """
    age = bin_age_obs              # только пространственная часть (num_bin,)
    err = bin_age_err_obs
    with np.errstate(divide="ignore", invalid="ignore"):
        cv        = np.where(age > 0, err / age, np.nan)
        sigma_log = np.sqrt(np.log1p(cv ** 2))
        mu_log    = np.where(age > 0, np.log(age) - 0.5 * sigma_log ** 2, np.nan)
    log_samples = np.random.normal(loc=mu_log, scale=sigma_log, size=(n_sim,) + age.shape)
    return np.exp(log_samples)   # (n_sim, num_bin)


def _mc_perturb_met_vec(bin_met_obs, bin_met_err_obs, n_sim):
    """
    Возмущение вектора металличности (нормальное в log-пространстве, сдвига нет).
    Возвращает массив (n_sim, num_bin).
    """
    return np.random.normal(
        loc=bin_met_obs,
        scale=bin_met_err_obs,
        size=(n_sim,) + bin_met_obs.shape,
    )   # (n_sim, num_bin)


def run_mc_stability(weight_matrix, n_sim=50):
    """
    Монте-Карло оценка устойчивости пространственных карт возраста и металличности.

    Алгоритм:
      1. Для каждой симуляции возмущаем наблюдаемые бин-векторы согласно ошибкам.
      2. Дополняем вектор нулями регуляризации (tail из num_bin-3 нулей).
      3. Решаем lsq_linear с теми же bounds и weight_matrix.
      4. Строим пространственную карту build_model_chem_image.
      5. Вычисляем mean и std по всем симуляциям.

    Возвращает:
        age_mean, age_std, met_mean, met_std — каждый (scale, scale).
    """
    print(f"MC-устойчивость: {n_sim} симуляций...")

    # Возмущённые бин-векторы (только пространственная часть)
    age_sims_vec = _mc_perturb_age_vec(bin_age[:num_bin], bin_age_err, n_sim)  # (n_sim, num_bin)
    met_sims_vec = _mc_perturb_met_vec(bin_met[:num_bin], bin_met_err, n_sim)

    age_images = []
    met_images = []

    for i in tqdm(range(n_sim)):
        # Вектор с регуляризацией: первые num_bin — возмущённые данные, хвост — нули
        b_age = np.zeros(num_bin * 2 - 3)
        b_met = np.zeros(num_bin * 2 - 3)
        b_age[:num_bin] = age_sims_vec[i]
        b_met[:num_bin] = met_sims_vec[i]

        phi_age_i = scipy.optimize.lsq_linear(
            weight_matrix, b_age,
            bounds=(age_min, age_max), method="bvls",
            tol=1e-10, max_iter=None, verbose=0,
        )
        phi_met_i = scipy.optimize.lsq_linear(
            weight_matrix, b_met,
            bounds=(met_min, met_max), method="bvls",
            tol=1e-10, max_iter=None, verbose=0,
        )

        age_images.append(build_model_chem_image(phi_age_i.x, weight_matrix))
        met_images.append(build_model_chem_image(phi_met_i.x, weight_matrix))

    age_stack = np.array(age_images)   # (n_sim, scale, scale)
    met_stack = np.array(met_images)

    return (
        np.nanmean(age_stack, axis=0), np.nanstd(age_stack, axis=0),
        np.nanmean(met_stack, axis=0), np.nanstd(met_stack, axis=0),
    )




def compute_r_covered(weight_matrix):
    """
    Возвращает максимальный радиус наблюдательного покрытия (в угловых секундах).

    Берётся напрямую из bin_scheme — координаты центров наблюдательных бинов
    уже хранятся в угловых секундах, поэтому максимальный радиус от центра даёт
    реальную границу шестиугольника IFU, независимо от орбитальной шкалы R_max.

    Для радиальных профилей: 21 бин равномерно укладывается в [0, r_covered],
    что гарантирует, что ось X профиля совпадает с видимой областью на картах.
    """
    # Координаты бинов в arcsec (столбцы 0 и 1 bin_scheme)
    r_bins = np.sqrt(bin_scheme[:, 0]**2 + bin_scheme[:, 1]**2)
    r_covered = np.max(r_bins)

    print(f"  Покрытие наблюдений (из bin_scheme): {r_covered:.2f}″"
          f"  (R_max орбит = {R_max:.2f}″)")
    return r_covered

def _compute_radial_profiles(model_age, model_met, weight_matrix, r_max=None):
    """
    Вычисляет радиальные профили возраста и металличности
    для всех кинематических подмножеств из LZ_REGIONS + All.

    Вес каждого динамического компонента — сумма физических весов орбит
    weights[model_index][orbi] по всем орбитам ячейки (idx_r, idx_lz).
    Это корректный физический вес (звёздная масса), не зависящий от нормировки
    матрицы weight_matrix. Компоненты с model_age = 0 (граница lsq_linear)
    получают правильно малый вес, если их орбитальная масса мала.

    Возвращает dict:
        { label: (age_profile_21, met_profile_21) }
    где профиль — массив (21,), NaN там где нет данных.
    """
    # Физические веса компонент: сумма weights по орбитам каждой ячейки
    comp_mass = np.zeros(len(dyn_comps_data))
    for comp_i, (idx_r, idx_lz) in enumerate(dyn_comps_data):
        orb_list = sorted_orbs[idx_r][idx_lz]
        if orb_list:
            comp_mass[comp_i] = np.sum(weights[model_index][orb_list])

    idx_r_all  = dyn_comps_data[:, 0]
    idx_lz_all = dyn_comps_data[:, 1]

    subsets = [(COMP_LABELS_RU.get("All", "Все компоненты"), np.ones(len(dyn_comps_data), dtype=bool))]
    for lz_start, lz_end, tag, _ in LZ_REGIONS:
        mask = (idx_lz_all >= lz_start) & (idx_lz_all < lz_end)
        subsets.append((COMP_LABELS_RU.get(tag, tag.capitalize()), mask))

    result = {}
    for label, lz_mask in subsets:
        age_prof = np.full(21, np.nan)
        met_prof = np.full(21, np.nan)
        for ir in range(21):
            sel = (idx_r_all == ir) & lz_mask
            if not np.any(sel):
                continue
            w = comp_mass[sel]
            w_sum = np.sum(w)
            if w_sum == 0:
                continue
            age_prof[ir] = np.sum(w * model_age[sel]) / w_sum
            met_prof[ir] = np.sum(w * model_met[sel]) / w_sum
        result[label] = (age_prof, met_prof)

    return result


def run_mc_radial_profiles(weight_matrix, n_sim=50, r_max=None):
    """
    Монте-Карло оценка неопределённости радиальных профилей.

    Для каждой симуляции:
      1. Возмущаем bin_age / bin_met согласно ошибкам.
      2. Решаем lsq_linear.
      3. Вычисляем радиальные профили для всех компонент.

    Возвращает dict:
        { label: (age_mean_21, age_std_21, met_mean_21, met_std_21) }
    """
    print(f"MC радиальные профили: {n_sim} симуляций...")

    age_sims_vec = _mc_perturb_age_vec(bin_age[:num_bin], bin_age_err, n_sim)
    met_sims_vec = _mc_perturb_met_vec(bin_met[:num_bin], bin_met_err, n_sim)

    # Накапливаем профили: {label: [список массивов (21,) по симуляциям]}
    age_accum = {}
    met_accum = {}

    for i in tqdm(range(n_sim)):
        b_age = np.zeros(num_bin * 2 - 3)
        b_met = np.zeros(num_bin * 2 - 3)
        b_age[:num_bin] = age_sims_vec[i]
        b_met[:num_bin] = met_sims_vec[i]

        phi_age_i = scipy.optimize.lsq_linear(
            weight_matrix, b_age,
            bounds=(age_min, age_max), method="bvls",
            tol=1e-10, max_iter=None, verbose=0,
        )
        phi_met_i = scipy.optimize.lsq_linear(
            weight_matrix, b_met,
            bounds=(met_min, met_max), method="bvls",
            tol=1e-10, max_iter=None, verbose=0,
        )

        profiles = _compute_radial_profiles(phi_age_i.x, phi_met_i.x, weight_matrix, r_max=r_max)
        for label, (ap, mp) in profiles.items():
            age_accum.setdefault(label, []).append(ap)
            met_accum.setdefault(label, []).append(mp)

    result = {}
    for label in age_accum:
        age_stack = np.array(age_accum[label])   # (n_sim, 21)
        met_stack = np.array(met_accum[label])
        result[label] = (
            np.nanmean(age_stack, axis=0), np.nanstd(age_stack, axis=0),
            np.nanmean(met_stack, axis=0), np.nanstd(met_stack, axis=0),
        )
    return result


def make_mc_radial_profiles_plot(
    model_age, model_met, weight_matrix,
    mc_profiles,
    r_max=None,
    filename="radial_profiles_mc.pdf",
):
    """
    Радиальные профили возраста и металличности с погрешностями из MC.

    Компоновка: N строк × 2 колонки, где каждая строка — один кинематический
    компонент (All / Corotating / Spherical / Counterrotating).
    Левая колонка — возраст, правая — металличность.

    В каждой ячейке:
      - Сплошная линия — номинальный профиль из основного решения.
      - Полупрозрачная полоса — ±1σ из MC (fill_between, step="mid").

    Оси X одинаковы для всех строк (sharex), подпись только снизу.
    Оси Y одинаковы в пределах каждой колонки (sharey="col").
    Название компонента — справа от строки через twinx.
    """
    _r_scale     = r_max if r_max is not None else R_max
    comp_nominal = _compute_radial_profiles(model_age, model_met, weight_matrix, r_max=_r_scale)
    r_arcsec     = np.arange(21) / 20.0 * _r_scale

    # Порядок строк: All первым, затем компоненты из LZ_REGIONS
    ordered_labels = (
        [COMP_LABELS_RU.get("All", "Все компоненты")] +
        [COMP_LABELS_RU.get(tag, tag.capitalize()) for _, _, tag, _ in LZ_REGIONS]
    )

    colors_map = {
        COMP_LABELS_RU.get("All", "Все компоненты"):             "black",
        COMP_LABELS_RU.get("corotating", "Совращающийся"):            "tab:blue",
        COMP_LABELS_RU.get("spherical",  "Сферический"):         "tab:green",
        COMP_LABELS_RU.get("counterrotating", "Противовращающийся"):   "tab:red",
    }

    N = len(ordered_labels)
    fig, axes = plt.subplots(
        nrows=N, ncols=2,
        figsize=(12, 3.5 * N),
        sharex=True,
        sharey="col",
        gridspec_kw={"hspace": 0.06, "wspace": 0.25},
    )

    for row, label in enumerate(ordered_labels):
        color = colors_map.get(label, "grey")
        ax_age = axes[row, 0]
        ax_met = axes[row, 1]

        age_nom, met_nom = comp_nominal[label]
        age_mean, age_std, met_mean, met_std = mc_profiles[label]

        for ax, nom, mean, std, ylabel in [
            (ax_age, age_nom, age_mean, age_std, "Возраст (млрд лет)"),
            (ax_met, met_nom, met_mean, met_std, "Металличность [Fe/H]"),
        ]:
            valid = ~np.isnan(nom)
            rv = r_arcsec[valid]

            # Номинальный профиль
            ax.step(rv, nom[valid],
                    where="mid", color=color, lw=2.2, ls="-")

            # ±1σ полоса вокруг номинала (std из MC — оценка неопределённости,
            # центр — номинальное решение, а не MC-среднее)
            ax.fill_between(rv,
                            (nom - std)[valid],
                            (nom + std)[valid],
                            step="mid", alpha=0.20, color=color, linewidth=0)

            ax.set_xlim(left=0)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=11)

            # Y-подпись только для первой строки
            if row == 0:
                ax.set_title(ylabel, fontsize=13, pad=6)

            # Y-метка слева
            ax.set_ylabel(ylabel, fontsize=11)

            # X-подпись только в последней строке
            if row == N - 1:
                ax.set_xlabel(r"$r$ (″)", fontsize=13)

        # Название компонента справа от строки
        ax_r = ax_met.twinx()
        ax_r.set_yticks([])
        ax_r.set_ylabel(label, fontsize=13, rotation=270, labelpad=20,
                        va="bottom", color=color)

    fig.savefig(filename, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Радиальные профили с MC-погрешностями сохранены: {filename}")

def make_mc_stability_plot(
    model_age_image, model_met_image,
    age_mean, age_std,
    met_mean, met_std,
    age_vmin, age_vmax,
    met_vmin, met_vmax,
    filename="mc_stability.pdf",
):
    """
    Сетка 2×4: две строки (возраст | металличность), четыре колонки:
      Исходная модель | MC-среднее | Разность (mean − model) | MC-std

    Единая палитра на всю фигуру: шкала и колормап одинаковы для всех
    ячеек одной строки. Границы вычисляются по nanmin/nanmax всех четырёх
    карт строки вместе.

    Колорбары:
      - Сверху над строкой «Возраст»     (ticks сверху)
      - Снизу под строкой «Металличность» (ticks снизу)
    """
    from matplotlib.gridspec import GridSpec
    import matplotlib.cm     as mplcm
    import matplotlib.colors as mplcolors

    age_diff = age_mean - model_age_image
    met_diff = met_mean - model_met_image

    # Единая шкала каждой строки: nanmin/nanmax по всем четырём картам
    def row_lim(*arrs):
        vals = np.concatenate([np.asarray(a).ravel() for a in arrs])
        return np.nanpercentile(vals, 2), np.nanpercentile(vals, 98)

    age_vmin_u, age_vmax_u = row_lim(model_age_image, age_mean, age_diff, age_std)
    met_vmin_u, met_vmax_u = row_lim(model_met_image, met_mean, met_diff, met_std)

    # rows_cfg: (карты × 4, cmap, vmin, vmax, подпись строки)
    rows_cfg = [
        (model_age_image, age_mean, age_diff, age_std,
         cmap_chem, age_vmin_u, age_vmax_u, "Возраст (млрд лет)"),
        (model_met_image, met_mean, met_diff, met_std,
         cmap_chem, met_vmin_u, met_vmax_u, "Металличность [Fe/H]"),
    ]

    COL_LABELS = ["Модель", "MC-среднее", "Разность (мод.−набл.)", "MC-разброс (σ)"]
    N_ROWS, N_COLS = 2, 4

    LEFT    = 0.10
    RIGHT   = 0.98
    CBAR_H  = 0.032
    GAP     = 0.005

    # Зоны (снизу вверх):
    #  BOTTOM_CB (колорбар металличности, снизу)
    #  карты
    #  TOP_CB    (колорбар возраста, сверху)
    BOTTOM_CB_BOT = 0.02
    BOTTOM_CB_TOP = BOTTOM_CB_BOT + CBAR_H
    MAPS_BOT      = BOTTOM_CB_TOP + GAP
    MAPS_TOP      = 0.94 - CBAR_H - GAP
    TOP_CB_BOT    = MAPS_TOP + GAP
    TOP_CB_TOP    = TOP_CB_BOT + CBAR_H

    fig = plt.figure(figsize=(16, 9))

    gs_maps    = GridSpec(N_ROWS, N_COLS, hspace=0, wspace=0,
                          left=LEFT, right=RIGHT,
                          bottom=MAPS_BOT, top=MAPS_TOP, figure=fig)
    gs_cb_top  = GridSpec(1, N_COLS, hspace=0, wspace=0,   # возраст — сверху
                          left=LEFT, right=RIGHT,
                          bottom=TOP_CB_BOT, top=TOP_CB_TOP, figure=fig)
    gs_cb_bot  = GridSpec(1, N_COLS, hspace=0, wspace=0,   # металличность — снизу
                          left=LEFT, right=RIGHT,
                          bottom=BOTTOM_CB_BOT, top=BOTTOM_CB_TOP, figure=fig)

    fmt = tcr.FuncFormatter(_arcsec_formatter)
    _px_zero     = scale / 2 - 0.5 * BIN_SCHEME_SCALE
    _half_arcsec = (scale / 2) * _arcsec_per_pixel
    _step_arcsec = min([1, 2, 5, 10, 15, 20, 25, 30],
                       key=lambda v: abs(v - _half_arcsec * 0.6))
    _step_px    = _step_arcsec * BIN_SCHEME_SCALE
    _sym3_ticks = [_px_zero - _step_px, _px_zero, _px_zero + _step_px]

    def _add_cbar(gs, col, cmap, vmin, vmax, label, ticks_pos):
        cax  = fig.add_subplot(gs[0, col])
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
        sm   = mplcm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb   = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(label, labelpad=4, fontsize=13)
        cb.ax.tick_params(labelsize=11)
        cax.xaxis.set_ticks_position(ticks_pos)
        cax.xaxis.set_label_position(ticks_pos)

    for ci, (_, _, _, _, cmap, vmin, vmax, row_label) in enumerate(rows_cfg):
        col_label = f"{row_label}  [{COL_LABELS[ci % N_COLS]}]" if False else COL_LABELS[ci]
    # Колорбар возраста — сверху (ticks top)
    for ci in range(N_COLS):
        _add_cbar(gs_cb_top, ci,
                  rows_cfg[0][4], rows_cfg[0][5], rows_cfg[0][6],
                  COL_LABELS[ci], "top")
    # Колорбар металличности — снизу (ticks bottom)
    for ci in range(N_COLS):
        _add_cbar(gs_cb_bot, ci,
                  rows_cfg[1][4], rows_cfg[1][5], rows_cfg[1][6],
                  COL_LABELS[ci], "bottom")

    # ── Карты ─────────────────────────────────────────────────────────────────
    for row, (m0, m1, m2, m3, cmap, vmin, vmax, row_label) in enumerate(rows_cfg):
        for col, data in enumerate([m0, m1, m2, m3]):
            ax = fig.add_subplot(gs_maps[row, col])
            ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax,
                      origin="lower", aspect="equal", interpolation="nearest")

            # X-ось: только нижняя строка
            if row == N_ROWS - 1:
                ax.xaxis.set_major_formatter(fmt)
                ax.xaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="x", labelsize=12, labelbottom=True)
                ax.set_xlabel("угл. сек.", fontsize=13)
            else:
                ax.tick_params(axis="x", labelbottom=False)
                ax.set_xlabel("")

            # Y-ось: только первая колонка
            if col == 0:
                ax.yaxis.set_major_formatter(fmt)
                ax.yaxis.set_major_locator(tcr.FixedLocator(_sym3_ticks))
                ax.tick_params(axis="y", labelsize=12, labelleft=True)
                ax.set_ylabel(row_label, fontsize=9)
            else:
                ax.tick_params(axis="y", labelleft=False)
                ax.set_ylabel("")

    fig.savefig(filename, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"MC-устойчивость сохранена: {filename}")

# ========== ОСНОВНОЙ КОД ==========

def main():
    # ── weight_matrix ─────────────────────────────────────────────────────────
    _wm_cache = _rload("weight_matrix")
    if _wm_cache is not None:
        weight_matrix    = _wm_cache["weight_matrix"]
        weight_map_binned = _wm_cache["weight_map_binned"]
    else:
        weight_matrix, weight_map_binned = create_weight_matrix_of_dyn_comps(
            lambda_reg=DEFAULT_LAMBDA_REG
        )
        _rsave("weight_matrix",
               weight_matrix=weight_matrix,
               weight_map_binned=weight_map_binned)

    # Массовые доли компонент и русские метки
    mass_fractions = compute_mass_fractions()
    COMP_LABELS_RU.update(make_comp_labels_ru(mass_fractions))
    print("Массовые доли кинематических компонент:")
    for tag, frac in mass_fractions.items():
        if tag != "All":
            print(f"  {COMP_LABELS_RU[tag]}: {frac*100:.1f}%")

    # R_COVERED: граница наблюдательного покрытия — вычисляется всегда быстро
    # из weight_matrix, поэтому не кэшируется отдельно и корректно работает
    # со старым кэшем weight_matrix.
    R_COVERED = compute_r_covered(weight_matrix)

    # ── lsq_linear: model_age / model_met ────────────────────────────────────
    print("=" * 30)
    _lsq_cache = _rload("lsq_solution")
    if _lsq_cache is not None:
        model_age = _lsq_cache["model_age"]
        model_met = _lsq_cache["model_met"]
        print("lsq_linear: загружено из кэша.")
    else:
        print("Решение системы уравнений с помощью метода наименьших квадратов.")
        print("Используется scipy.optimize.lsq_linear, method='bvls':")
        phi_age = scipy.optimize.lsq_linear(
            weight_matrix, bin_age,
            bounds=(age_min, age_max), method='bvls',
            tol=1e-30, lsq_solver=None, lsmr_tol=None,
            max_iter=None, verbose=1, lsmr_maxiter=None
        )
        phi_met = scipy.optimize.lsq_linear(
            weight_matrix, bin_met,
            bounds=(met_min, met_max), method='bvls',
            tol=1e-30, lsq_solver=None, lsmr_tol=None,
            max_iter=None, verbose=1, lsmr_maxiter=None
        )
        print("Результат lsq_linear для возрастов:"); print(phi_age)
        print("Результат lsq_linear для металличности:"); print(phi_met)
        model_age = phi_age.x
        model_met = phi_met.x
        _rsave("lsq_solution", model_age=model_age, model_met=model_met)

    # ── Пространственные карты химии ──────────────────────────────────────────
    _chem_cache = _rload("chem_images")
    if _chem_cache is not None:
        model_age_image = _chem_cache["model_age_image"]
        model_met_image = _chem_cache["model_met_image"]
        age_by_tag = {tag: _chem_cache[f"age_{tag}"] for _, _, tag, _ in LZ_REGIONS}
        met_by_tag = {tag: _chem_cache[f"met_{tag}"] for _, _, tag, _ in LZ_REGIONS}
    else:
        model_age_image = build_model_chem_image(model_age, weight_matrix)
        model_met_image = build_model_chem_image(model_met, weight_matrix)
        age_by_tag = {}
        met_by_tag = {}
        for lz_start, lz_end, tag, _ in LZ_REGIONS:
            age_by_tag[tag] = build_model_chem_image_lz_region(
                model_age, weight_matrix, lz_start, lz_end)
            met_by_tag[tag] = build_model_chem_image_lz_region(
                model_met, weight_matrix, lz_start, lz_end)
        _rsave("chem_images",
               model_age_image=model_age_image,
               model_met_image=model_met_image,
               **{f"age_{tag}": age_by_tag[tag] for _, _, tag, _ in LZ_REGIONS},
               **{f"met_{tag}": met_by_tag[tag] for _, _, tag, _ in LZ_REGIONS})

    # ── Динамические цветовые границы ─────────────────────────────────────────
    age_vmin, age_vmax = _joint_limits(model_age_image, age_map)
    met_vmin, met_vmax = _joint_limits(model_met_image, met_map)
    print("=" * 30)
    print(f"Динамические границы шкалы возраста:       [{age_vmin:.3f}, {age_vmax:.3f}] млрд лет")
    print(f"Динамические границы шкалы металличности:  [{met_vmin:.3f}, {met_vmax:.3f}] dex")

    # ── Кинематические карты ──────────────────────────────────────────────────
    _kin_cache = _rload("kinem_maps")
    if _kin_cache is not None:
        velocity_map = _kin_cache["velocity_map"]
        sigma_map    = _kin_cache["sigma_map"]
        kin_by_tag   = {tag: (_kin_cache[f"vel_{tag}"], _kin_cache[f"sig_{tag}"])
                        for _, _, tag, _ in LZ_REGIONS}
    else:
        kin_all      = DYN_COMP_LOSVD_MAP_INDICES(0, 21)
        velocity_map = kin_all[0]
        sigma_map    = kin_all[1]
        kin_by_tag   = {}
        for lz_start, lz_end, tag, _ in LZ_REGIONS:
            kin = DYN_COMP_LOSVD_MAP_INDICES(lz_start, lz_end)
            kin_by_tag[tag] = (kin[0], kin[1])
        _rsave("kinem_maps",
               velocity_map=velocity_map, sigma_map=sigma_map,
               **{f"vel_{tag}": kin_by_tag[tag][0] for _, _, tag, _ in LZ_REGIONS},
               **{f"sig_{tag}": kin_by_tag[tag][1] for _, _, tag, _ in LZ_REGIONS})

    print(np.shape(velocity_map), np.shape(model_age))
    print(velocity_map)
    print(model_age)

    # ── Отдельные PDF (прежнее поведение) ────────────────────────────────────
    make_kinem_plot(velocity_map, label=r"$V$ (км/с)",       filename=f"{GALAXY}_velmap_model_all.pdf",  transpose=True)
    make_kinem_plot(sigma_map,    label=r"$\sigma$ (км/с)", filename=f"{GALAXY}_sigmap_model_all.pdf",  transpose=True)
    make_kinem_plot(velocity_obs, label=r"$V_\mathrm{набл}$ (км/с)",       filename=f"{GALAXY}_velmap_obs.pdf")
    make_kinem_plot(sigma_obs,    label=r"$\sigma_\mathrm{набл}$ (км/с)", filename=f"{GALAXY}_sigmap_obs.pdf")

    for lz_start, lz_end, tag, lz_label in LZ_REGIONS:
        make_kinem_plot(kin_by_tag[tag][0], label=rf"$V$ (км/с),  {lz_label}",
                        filename=f"{GALAXY}_velmap_model_{tag}.pdf", transpose=True)
        make_kinem_plot(kin_by_tag[tag][1], label=rf"$\sigma$ (км/с),  {lz_label}",
                        filename=f"{GALAXY}_sigmap_model_{tag}.pdf", transpose=True)

    make_lambda_z_vs_radius_vs_chem_plot(model_age, dyn_comps_data, age_vmin, age_vmax, chem="age")
    make_lambda_z_vs_radius_vs_chem_plot(model_met, dyn_comps_data, met_vmin, met_vmax, chem="met")
    make_obs_chem_plot(age_map,         age_vmin, age_vmax, chem="age")
    make_obs_chem_plot(met_map,         met_vmin, met_vmax, chem="met")
    make_model_chem_plot(model_age_image, age_vmin, age_vmax, chem="age_all")
    make_model_chem_plot(model_met_image, met_vmin, met_vmax, chem="met_all")
    for lz_start, lz_end, tag, lz_label in LZ_REGIONS:
        make_model_chem_plot(age_by_tag[tag], age_vmin, age_vmax, chem=f"age_{tag}")
        make_model_chem_plot(met_by_tag[tag], met_vmin, met_vmax, chem=f"met_{tag}")

    make_model_dencity_plot(weight_map_binned)

    # ── Сетка по металличности ────────────────────────────────────────────────
    make_met_grid(
        weight_matrix=weight_matrix,
        weight_map_binned=weight_map_binned,
        model_met=model_met,
        met_vmin=met_vmin,
        met_vmax=met_vmax,
        n_subsets=4,
        filename=f"{GALAXY}_met_grid.pdf",
    )

    make_age_grid(
        weight_matrix=weight_matrix,
        weight_map_binned=weight_map_binned,
        model_age=model_age,
        age_vmin=age_vmin,
        age_vmax=age_vmax,
        n_subsets=4,
        filename=f"{GALAXY}_age_grid.pdf",
    )

    # ── Радиальные профили химии ─────────────────────────────────────────────
    make_radial_chem_profiles(
        model_age=model_age,
        model_met=model_met,
        weight_matrix=weight_matrix,
        r_max=R_COVERED,
        filename=f"{GALAXY}_radial_chem_profiles.pdf",
    )

    make_radial_chem_profiles_stacked(
        model_age=model_age,
        model_met=model_met,
        weight_matrix=weight_matrix,
        r_max=R_COVERED,
        filename=f"{GALAXY}_radial_chem_profiles_stacked.pdf",
    )

    # ── MC-устойчивость ──────────────────────────────────────────────────────
    N_MC = 50

    # MC карты устойчивости
    _mc_stab_cache = _rload("mc_stability")
    if _mc_stab_cache is not None:
        age_mc_mean = _mc_stab_cache["age_mc_mean"]
        age_mc_std  = _mc_stab_cache["age_mc_std"]
        met_mc_mean = _mc_stab_cache["met_mc_mean"]
        met_mc_std  = _mc_stab_cache["met_mc_std"]
    else:
        age_mc_mean, age_mc_std, met_mc_mean, met_mc_std = run_mc_stability(
            weight_matrix, n_sim=N_MC,
        )
        _rsave("mc_stability",
               age_mc_mean=age_mc_mean, age_mc_std=age_mc_std,
               met_mc_mean=met_mc_mean, met_mc_std=met_mc_std)

    make_mc_stability_plot(
        model_age_image=model_age_image, model_met_image=model_met_image,
        age_mean=age_mc_mean, age_std=age_mc_std,
        met_mean=met_mc_mean, met_std=met_mc_std,
        age_vmin=age_vmin, age_vmax=age_vmax,
        met_vmin=met_vmin, met_vmax=met_vmax,
        filename=f"{GALAXY}_mc_stability.pdf",
    )

    # MC радиальные профили
    # ВНИМАНИЕ: если поменялась логика весов comp_weights, удалите mc_radial_profiles.npz
    _mc_rad_cache = _rload("mc_radial_profiles")
    if _mc_rad_cache is not None:
        # Таблица перевода старых (английских) ключей в текущие (русские).
        # Нужна для совместимости с кэшем, сохранённым до смены меток.
        _key_remap = {
            "All":             COMP_LABELS_RU.get("All",             "Все компоненты"),
            "Corotating":      COMP_LABELS_RU.get("corotating",      "Совращающийся"),
            "Spherical":       COMP_LABELS_RU.get("spherical",       "Сферический"),
            "Counterrotating": COMP_LABELS_RU.get("counterrotating", "Противовращающийся"),
        }
        mc_radial = {}
        labels = [s.decode() if isinstance(s, bytes) else s
                  for s in _mc_rad_cache["labels"]]
        for label in labels:
            new_key = _key_remap.get(label, label)   # старый → новый; неизвестные — как есть
            mc_radial[new_key] = (
                _mc_rad_cache[f"{label}_age_mean"],
                _mc_rad_cache[f"{label}_age_std"],
                _mc_rad_cache[f"{label}_met_mean"],
                _mc_rad_cache[f"{label}_met_std"],
            )
    else:
        mc_radial = run_mc_radial_profiles(weight_matrix, n_sim=N_MC, r_max=R_COVERED)
        _rsave("mc_radial_profiles",
               labels=np.array(list(mc_radial.keys())),
               **{f"{lbl}_age_mean": v[0] for lbl, v in mc_radial.items()},
               **{f"{lbl}_age_std":  v[1] for lbl, v in mc_radial.items()},
               **{f"{lbl}_met_mean": v[2] for lbl, v in mc_radial.items()},
               **{f"{lbl}_met_std":  v[3] for lbl, v in mc_radial.items()})

    make_mc_radial_profiles_plot(
        model_age=model_age, model_met=model_met,
        weight_matrix=weight_matrix,
        mc_profiles=mc_radial,
        r_max=R_COVERED,
        filename=f"{GALAXY}_radial_profiles_mc.pdf",
    )

    # ── Карты невязок ────────────────────────────────────────────────────────
    make_residuals_grid(
        model_age_image=model_age_image, age_map=age_map,
        model_met_image=model_met_image, met_map=met_map,
        velocity_map=velocity_map,       velocity_obs=velocity_obs,
        sigma_map=sigma_map,             sigma_obs=sigma_obs,
        age_vmin=age_vmin, age_vmax=age_vmax,
        met_vmin=met_vmin, met_vmax=met_vmax,
        filename=f"{GALAXY}_residuals_grid.pdf",
    )

    # ── Сводная таблица ───────────────────────────────────────────────────────
    # Порядок строк таблицы: obs | all | corotating | spherical | counterrotating
    TABLE_ROW_TAGS = ["corotating", "spherical", "counterrotating"]

    age_list = [age_map,         model_age_image] + [age_by_tag[t] for t in TABLE_ROW_TAGS]
    met_list = [met_map,         model_met_image] + [met_by_tag[t] for t in TABLE_ROW_TAGS]
    vel_list = [velocity_obs,    velocity_map]    + [kin_by_tag[t][0] for t in TABLE_ROW_TAGS]
    sig_list = [sigma_obs,       sigma_map]       + [kin_by_tag[t][1] for t in TABLE_ROW_TAGS]

    # obs — без транспонирования; всё остальное (модель) — с .T
    vel_tr = [False, True, True, True, True]
    sig_tr = [False, True, True, True, True]

    make_summary_grid(
        age_images=age_list, met_images=met_list,
        vel_maps=vel_list,   vel_transpose=vel_tr,
        sig_maps=sig_list,   sig_transpose=sig_tr,
        age_vmin=age_vmin,   age_vmax=age_vmax,
        met_vmin=met_vmin,   met_vmax=met_vmax,
        filename=f"{GALAXY}_summary_grid.pdf",
    )


# ========== ЗАПУСК ==========
if __name__ == "__main__":
    main()
