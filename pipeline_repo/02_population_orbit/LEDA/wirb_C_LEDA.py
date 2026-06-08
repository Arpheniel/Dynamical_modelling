import numpy as np
import matplotlib.pyplot as plt
from matplotlib import transforms, ticker as tcr, colors, patheffects as pe
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from scipy.ndimage import rotate as ndi_rotate
import agama
import os

# ========== КОНСТАНТЫ И НАСТРОЙКИ LEDA 2220522 ==========
scale = 32  # Размер картинки для LEDA (из neworbmap.py)
arcsec_to_kpc = 0.5296  # Из neworbmap.py
N = 20000        # финал C: 20000 орбит
M_total = 1.5e9  # светимость-нормировка подписи колорбара (×Υ = масса)
distance = 108643000  # Из run-скрипта LEDA: 108643 кпк
# Реальное наклонение галактики (из run-script, INCL = beta).
real_inclination = 27

# Изображаемое наклонение для каждого столбца плотности (в градусах).
# Можно поставить любое значение от 0 (face-on, вид с полюса) до 90
# (edge-on, ребро). Значение совпадает с тем, что будет на подписи.
displayed_inclination_col0 = real_inclination   # col 0: вид как с Земли
displayed_inclination_col1 = 90                 # col 1: edge-on

# incl_wirb — угол поворота rotation_matrix_x в orbit_map. Agama хранит
# orbit-positions в galactic frame (диск в плоскости xy), поэтому
# incl_wirb напрямую равно желаемому изображаемому наклонению.
incl_col0 = displayed_inclination_col0
incl_col1 = displayed_inclination_col1

# Для совместимости со старым кодом: incl ниже — для col 1.
# col 0 теперь использует incl_col0 явно.
incl = incl_col1
model_index = 5  # финал C (приор v=125): Υ=4.781, индекс 5 из 9
gamma = 0  # ФИКС: лишний self.tr-поворот убран. PA задаётся через agama_gamma
           # (ниже, со знаком −). Иначе поворот двигал только кинематику → плотность ⟂ скорость.

# ========== ПАРАМЕТРЫ AGAMA-ПРОЕКЦИИ ==========
# Углы Эйлера, использованные при построении модели в run_forstand_grid_*.py.
# Должны точно совпадать с тем, что было передано в Target(type='LOSVD',
# alpha=..., beta=..., gamma=...) — иначе плотность не совпадёт по ориентации
# с кинематикой и наблюдениями.
agama_alpha = 0.0    # ALPHA из run-скрипта (0 для LEDA)
agama_beta  = 27.0   # INCL из run-скрипта (27° для LEDA, наклонение)
agama_gamma = -40.0  # ФИКС: знак инвертирован (agama-углы разнознаковы, как beta в edge-on).
                     # +40 давал плотность ⟂ кинематике (.T=130° vs скорость 50°); −40 → .T=50° (совпало).
                     # gamma1 модели = 40°, PA на небе.
                     # Должно совпадать с тем, что было передано в
                     # Target(type='LOSVD', gamma=...) при построении модели.

# Угол поворота МОДЕЛЬНОЙ ПЛОТНОСТИ относительно горизонтали.
#
# С момента перевода проекции на (X, Y) (плоскость неба) поворот gamma1
# уже зашит в координаты orbit-positions (Agama сохраняет их в observer
# frame). Поэтому density_PA = 0 даёт ФИЗИЧЕСКИ КОРРЕКТНЫЙ результат —
# никакой постпроцессинговой поправки не нужно.
#
# Параметр оставлен на случай, если понадобится подкрутить вручную
# (например, при отладке или если Agama хранит orbits в другой системе).
density_PA = 0.0

# Преобразования модельной плотности для согласования с кинематикой.
# С (X, Y)-проекцией плотность строится прямо в плоскости неба (как и
# кинематика через bin_scheme), поэтому никаких преобразований не нужно.
density_transpose = True  # компенсация .T кинематики через flip_diagonal=True

# Доп. flip по X (если после transpose осталось зеркальное смещение).
# По умолчанию False; включите, если плотность отражена по X относительно
# кинематики и наблюдений.
density_flip_x = False

# Доп. flip по Y (аналогично).
density_flip_y = False

# ========== РАЗБИЕНИЕ ОРБИТ НА ТРИ КОМПОНЕНТЫ ==========
# Границы по индексам бинов λ_z (схема А, 0..20). Меняйте значения, чтобы
# попробовать разные варианты разбиения БЕЗ пересчёта модели — нужен только
# повторный запуск этого скрипта.
#
# По умолчанию симметричное разбиение около |λ_z| ≈ 0.35:
#   counter_rot:  индексы [0, 7)   ↔  λ_z ∈ [-1.00, -0.35]
#   spheroid:     индексы [7, 14)  ↔  λ_z ∈ (-0.35, +0.35)
#   co_rot:       индексы [14, 21) ↔  λ_z ∈ [+0.35, +1.00]
#
# Примеры альтернативных разбиений (со ссылкой на λ_z-границы):
#   ±0.25 (широкие диски):  counter_rot_end_idx=8,  co_rot_start_idx=13
#   ±0.35 (по умолчанию):   counter_rot_end_idx=7,  co_rot_start_idx=14
#   ±0.45 (узкие диски):    counter_rot_end_idx=6,  co_rot_start_idx=15
#   ±0.55 (очень узкие):    counter_rot_end_idx=5,  co_rot_start_idx=16
#   асимметрично (диск λ_z≥0.15, контр λ_z≤-0.45):
#                           counter_rot_end_idx=6,  co_rot_start_idx=12
#
# Симметрия около нуля ⇔ counter_rot_end_idx + co_rot_start_idx == 21.
#
# Для подписей на картинке («1.0 > λ_z > 0.35» и т.п.) границы вычисляются
# автоматически через bin_index_to_lambda_z_range().
counter_rot_end_idx = 7    # контр-вращение: range(0, counter_rot_end_idx)
co_rot_start_idx    = 14   # со-вращение:    range(co_rot_start_idx, 21)
# Сферическая компонента — всё между: range(counter_rot_end_idx, co_rot_start_idx)

Ph_min=0.0000000000001
Phot_min=0.000000001
# Цветовые карты (оставляем из первого кода)
cmap = plt.cm.jet
cmap.set_bad(color='White')
# ИЗМЕНЕНО: Черный фон для плохих значений на всех картах плотности
cmap_2 = plt.cm.gist_heat.copy()
cmap_2.set_bad(color='black')
cmap_2.set_under(color='black')
cmap.set_under(color="white")

# Черный фон для повернутой плотности
cmap_density_rotated = plt.cm.gist_heat.copy()
cmap_density_rotated.set_bad(color='black')
cmap_density_rotated.set_under(color='black')

# Лимиты скорости и дисперсии вычисляются автоматически по данным —
# см. вызов compute_global_kinematic_limits() ниже. Симметрично для скорости,
# расширенно по фактическим значениям для дисперсии.
Vel_max = Vel_min = None
Sigma_max = Sigma_min = None
Ph_max = 0.024  # Из второго кода

# Шрифты
fontdict = {'fontsize': 6}
font_italic = FontProperties(family='italic')

# Путь для кэширования карт плотности
CACHE_DIR_ORIGINAL = "density_cache_original_LEDA_agama_21r"  # проекция (X, Y)
CACHE_DIR_ROTATED = "density_cache_rotated_LEDA_agama_21r"   # проекция (X, Y)
for cache_dir in [CACHE_DIR_ORIGINAL, CACHE_DIR_ROTATED]:
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

# ========== ЗАГРУЗКА ДАННЫХ ==========
print("Загрузка данных...")
# ИЗМЕНЕНО: Имя файла из второго кода
Filename = "M1e+07_O0_Rh30.6_Vh125_i27_a0_N20000_R1.00_GH_DensityCylindricalLinear.npz"
archive = np.load(Filename, allow_pickle=True, encoding='latin1')

weights = archive["weights"]
ML = archive["Upsilon"][model_index]
Rmean_list = archive["Rmean"]

lambda_z_list = archive["MOD_Lambda_z"]
DYN_COMP_LOSVD = archive["DYN_COMP_LOSVD"][model_index]

# Размерность LOSVD (число апертур × число b-сплайнов на сетке скоростей).
# Берётся автоматически из первой непустой ячейки DYN_COMP_LOSVD —
# чтобы скрипт работал для любой галактики (у PGC 676 апертур, у LEDA 263).
LOSVD_SHAPE = None
for _i in range(21):
    for _j in range(21):
        _cell = DYN_COMP_LOSVD[_i][_j]
        if isinstance(_cell, np.ndarray) and _cell.size > 0:
            LOSVD_SHAPE = _cell.shape
            break
    if LOSVD_SHAPE is not None:
        break
if LOSVD_SHAPE is None:
    raise RuntimeError("DYN_COMP_LOSVD пустой — нечего обрабатывать")
print(f"  Размерность LOSVD: {LOSVD_SHAPE} (апертур × b-сплайнов)")

# ИЗМЕНЕНО: Файлы спектра, фотометрии и bin_scheme из второго кода
spectr = fits.open("results_8254-1902_vorb020_md19_ad-1_nmom4.fits")
phot = fits.open("mosaic-00157902-LEDA2220522-z-CCD3-image.fits")
bin_scheme = np.loadtxt("bins_LEDA_2220522_Damir's.txt")

# ========== ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ==========
def tickers_X_formatter(x, pos):
    """Форматтер для подписей в угловых секундах (с шагом 3) — для LEDA 2220522"""
    value = int((x - scale/2)/2)
    if value % 3 == 0 and abs(value) <= 6:
        return f'{value}'
    return ''

def tickers_Y_formatter(y, pos):
    """Форматтер для подписей в угловых секундах по Y (с шагом 3) — для LEDA 2220522"""
    value = int((y - scale/2)/2)
    if value % 3 == 0 and abs(value) <= 6:
        return f'{value}'
    return ''

def make_kpc_scale(scale_size, offset=0):
    """Создание шкалы в кпк из второго кода"""
    kpc_scale = np.zeros((1, int((scale_size/2)*arcsec_to_kpc+1)))[0]
    for kpc in range(0, len(kpc_scale)):
        kpc_scale[kpc] = kpc/arcsec_to_kpc * 2 + 3.45 + offset
    return kpc_scale

def tickers_X_formatter_kpc(x, pos):
    """Форматтер для подписей в килопарсеках"""
    return f'{abs((x - scale/2)/2 * arcsec_to_kpc):.0f}'

def DYN_COMP_LOSVD_MAP_INDICES(lz_start_idx, lz_end_idx):
    """
    Создание карт кинематики для заданного диапазона индексов λ_z
    lz_start_idx: начальный индекс (включительно)
    lz_end_idx: конечный индекс (не включая)
    """
    # ИЗМЕНЕНО: Размеры из второго кода (676 вместо 263)
    LOSVD_plt = np.zeros(LOSVD_SHAPE)
    
    # ИЗМЕНЕНО: Суммируем по заданным индексам
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
        # ИЗМЕНЕНО: Индексация bin_scheme из второго кода
        for bin_i in range(len(bin_scheme)):
            x = int(bin_scheme[bin_i][0] + scale/2)
            y = int(bin_scheme[bin_i][1] + scale/2)
            if 0 <= x < scale and 0 <= y < scale:
                KINEM_MAP[y][x] = GH_moments[int(bin_scheme[bin_i][2])][gh_moment]
        LOSVD_MAP.append(KINEM_MAP)
    
    return LOSVD_MAP

def lambda_z_to_index(lambda_z):
    """Преобразование λ_z в индекс (схема А: 21 бин с центрами 0..20).

    Используется int(round(...)) — согласовано с биннингом в schwarzlib
    (DYN_COMP_LOSVD и MOD_weights_mult_matrices). Бины:
      idx=0:  λ_z ∈ [-1.000, -0.950)   ширина 0.05
      idx=1:  λ_z ∈ [-0.950, -0.850)   ширина 0.10
      ...
      idx=10: λ_z ∈ [-0.050, +0.050)   ширина 0.10  (центральный, через 0)
      ...
      idx=19: λ_z ∈ [+0.850, +0.950)   ширина 0.10
      idx=20: λ_z ∈ [+0.950, +1.000]   ширина 0.05

    Внимание: крайние бины (0 и 20) физически вдвое уже центральных.
    """
    return int(round((1 + lambda_z) * 10))


def bin_index_to_lambda_z_range(start_idx, end_idx):
    """Возвращает (lo, hi) в λ_z для диапазона индексов [start_idx, end_idx).

    Используется для автоматических подписей на картинке: какому диапазону
    λ_z соответствует диапазон индексов компоненты.
    """
    lo = max((start_idx - 0.5) / 10 - 1, -1.0)
    hi = min((end_idx - 1 + 0.5) / 10 - 1, 1.0)
    return lo, hi


# Диапазоны индексов трёх компонент. Вычисляются из counter_rot_end_idx
# и co_rot_start_idx (см. начало файла) — поэтому, чтобы попробовать
# другое разбиение, достаточно поменять там два числа.
COMPONENT_LZ_RANGES = [
    (0, counter_rot_end_idx),               # контр-вращающийся диск
    (counter_rot_end_idx, co_rot_start_idx),# сферическая компонента
    (co_rot_start_idx, 21),                 # со-вращающийся диск
]


def compute_global_kinematic_limits(component_ranges, padding=1.05, round_to=10):
    """Вычисляет общие лимиты для карт скорости и дисперсии так, чтобы все
    значения, реально появляющиеся на картах всех компонент, попадали в шкалу.

    Скорость симметризуется: Vel_min = -Vel_lim, Vel_max = +Vel_lim.
    Дисперсия — расширяется в обе стороны по факту.
    Округление до значений, кратных round_to (по умолчанию 10), для красивых
    чисел на колорбаре.

    Параметры:
        component_ranges: список пар (lz_start, lz_end) для каждой компоненты
        padding: множитель запаса (1.05 = +5%, чтобы данные не упирались в край)
        round_to: округлять лимиты до ближайшего значения, кратного этому числу

    Возвращает (Vel_lim, Sigma_min_val, Sigma_max_val).
    """
    # Берём только апертуры, которые реально отрисовываются (из bin_scheme),
    # чтобы пустые/неиспользуемые строки GH_moments не растягивали шкалу.
    used_aper = np.unique(bin_scheme[:, 2].astype(int))

    all_v, all_sigma = [], []
    for lz_start, lz_end in component_ranges:
        LOSVD_plt = np.zeros(LOSVD_SHAPE)
        for r in range(21):
            for l_z in range(lz_start, lz_end):
                LOSVD_plt += DYN_COMP_LOSVD[l_z][r]
        GH_moments = agama.ghMoments(
            matrix=LOSVD_plt * ML**-0.5,
            gridv=np.linspace(-250, 250, 46) * ML**0.5,
            degree=2, ghorder=6
        )[:, (1, 2, 6, 7, 8, 9)]

        v     = GH_moments[used_aper, 0]
        sigma = GH_moments[used_aper, 1]
        all_v.append(v[np.isfinite(v)])
        all_sigma.append(sigma[np.isfinite(sigma)])

    all_v     = np.concatenate(all_v)
    all_sigma = np.concatenate(all_sigma)

    # Симметричная шкала по v: ±max(|v|) * padding, округлённое вверх до round_to
    v_extreme = np.max(np.abs(all_v)) * padding
    Vel_lim = int(np.ceil(v_extreme / round_to) * round_to)

    # Расширенная шкала по σ
    Sigma_min_val = int(np.floor(all_sigma.min() / padding / round_to) * round_to)
    Sigma_max_val = int(np.ceil (all_sigma.max() * padding / round_to) * round_to)
    Sigma_min_val = max(Sigma_min_val, 0)  # σ — положительная

    print(f"  Автолимиты кинематики:  V ∈ [-{Vel_lim}, +{Vel_lim}] км/с,  "
          f"σ ∈ [{Sigma_min_val}, {Sigma_max_val}] км/с")
    return Vel_lim, Sigma_min_val, Sigma_max_val


print("Расчёт кинематических лимитов...")
_Vel_lim, Sigma_min, Sigma_max = compute_global_kinematic_limits(COMPONENT_LZ_RANGES)
Vel_min, Vel_max = -_Vel_lim, +_Vel_lim


# ========== ОРБИТАЛЬНОЕ МОДЕЛИРОВАНИЕ С КЭШИРОВАНИЕМ ==========
print("Сортировка орбит...")
# ИЗМЕНЕНО: 21x21 (радиус x λ_z) - возвращаем оригинальную размерность
sorted_orbs = [[[] for _ in range(21)] for _ in range(21)]

for orbi in range(N):
    Rmean = Rmean_list[orbi]
    lambda_z = lambda_z_list[orbi]
    # Схема А: int(round(...)) — те же индексы, что в schwarzlib
    # (раньше здесь стояло int(...) без round, что давало рассогласование
    # карт плотности и кинематики на пограничных орбитах).
    idx_r = int(round((Rmean / (np.mean(Rmean_list) * 3)) * 20))
    idx_lz = int(round((1 + lambda_z) * 10))

    if 0 <= idx_r <= 20 and 0 <= idx_lz <= 20:
        sorted_orbs[idx_r][idx_lz].append(orbi)

def rotation_matrix_x(a):
    return np.array([[1, 0, 0],
                     [0, np.cos(a), -np.sin(a)],
                     [0, np.sin(a), np.cos(a)]])

def rotation_matrix_2d(a):
    """Матрица поворота в 2D плоскости"""
    c = np.cos(np.radians(a))
    s = np.sin(np.radians(a))
    return np.array([[c, -s],
                     [s, c]])


def rotate_density_nan_safe(image, angle_deg):
    """Поворот 2D-карты плотности на угол angle_deg вокруг центра, корректно
    обрабатывая NaN (без интерполяционного «расплывания»). Размер выходного
    массива совпадает с входным (нет смещения масштаба), при этом
    исходная картинка предварительно дополняется паддингом, чтобы поворот
    не обрезал угловые пиксели.

    Алгоритм:
      1) Заменяем NaN на 0 для самого изображения, готовим маску.
      2) Расширяем массив паддингом (sqrt(2) от размера) — чтобы при повороте
         ни один пиксель не вышел за пределы.
      3) Поворачиваем bilinear-интерполяцией (плавно).
      4) Поворачиваем маску nearest-neighbour, чтобы NaN не размывалась.
      5) Центрально обрезаем результат до исходного размера.
      6) Восстанавливаем NaN там, где маска её требует.
    """
    if angle_deg == 0:
        return image
    nan_mask = ~np.isfinite(image)
    img_filled = np.where(nan_mask, 0.0, image)

    # Паддинг, чтобы не было обрезки при повороте. sqrt(2)/2 ≈ 0.707 от
    # размера достаточно для любого угла; берём с запасом 1 пиксель.
    h, w = image.shape
    pad = int(np.ceil(max(h, w) * (np.sqrt(2) - 1) / 2)) + 1
    img_padded = np.pad(img_filled, pad, mode='constant', constant_values=0.0)
    mask_padded = np.pad(nan_mask.astype(np.float32), pad,
                         mode='constant', constant_values=1.0)

    rotated = ndi_rotate(img_padded, angle_deg, reshape=False,
                         order=1, cval=0.0)
    rotated_mask = ndi_rotate(mask_padded, angle_deg, reshape=False,
                              order=0, cval=1.0)

    # Central crop обратно до исходного размера.
    rotated = rotated[pad:pad + h, pad:pad + w]
    rotated_mask = rotated_mask[pad:pad + h, pad:pad + w]

    rotated[rotated_mask > 0.5] = np.nan
    return rotated


def detect_PA_from_velocity(v_map, threshold_frac=0.1):
    """Грубая оценка позиционного угла большой оси из карты лучевой скорости.

    Метод: для каждого угла θ ∈ [-90°, +90°) проецируем все пиксели на ось
    с этим углом и фитируем линейную V(s) = a·s + b. PA — это угол с
    максимальным |slope|·|r|, где r — коэффициент корреляции.

    threshold_frac: учитываются только пиксели с |V| > threshold_frac · max(|V|),
    чтобы шум вокруг нуля не размывал оценку.

    Возвращает PA в градусах в матплотлибовской конвенции (CCW от +x).
    Если данных недостаточно — возвращает None.
    """
    ny, nx = v_map.shape
    yy, xx = np.indices(v_map.shape).astype(float)
    xx -= nx / 2.0
    yy -= ny / 2.0

    valid = np.isfinite(v_map)
    if valid.sum() < 50:
        return None
    v_max = np.nanmax(np.abs(v_map))
    valid &= (np.abs(v_map) > threshold_frac * v_max)
    if valid.sum() < 50:
        return None

    x_arr = xx[valid]
    y_arr = yy[valid]
    v_arr = v_map[valid]
    v_centered = v_arr - v_arr.mean()

    best_theta_deg, best_score = 0.0, -np.inf
    for theta_deg in np.arange(-90, 90, 1):
        theta = np.radians(theta_deg)
        proj = x_arr * np.cos(theta) + y_arr * np.sin(theta)
        proj_centered = proj - proj.mean()
        denom = np.sum(proj_centered ** 2)
        if denom <= 0:
            continue
        slope = np.sum(proj_centered * v_centered) / denom
        # коэф. корреляции
        v_var = np.sum(v_centered ** 2)
        if v_var <= 0:
            continue
        r = (np.sum(proj_centered * v_centered)
             / np.sqrt(denom * v_var))
        score = abs(slope) * abs(r)
        if score > best_score:
            best_score = score
            best_theta_deg = theta_deg
    return float(best_theta_deg)


def get_cache_filename(lambda_z_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    """Генерация имени файла для кэша"""
    lz_min, lz_max = lambda_z_range
    if plane_rotation == 0:
        return os.path.join(cache_dir, f"density_lz_{lz_min}_{lz_max}_incl_{incl}.npy")
    else:
        return os.path.join(cache_dir, f"density_lz_{lz_min}_{lz_max}_incl_{incl}_plane_rot_{plane_rotation}.npy")

def save_density_map(density_map, lambda_z_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    """Сохранение карта плотности в файл"""
    filename = get_cache_filename(lambda_z_range, incl, plane_rotation, cache_dir)
    np.save(filename, density_map)
    print(f"  Карта сохранена в {filename}")
    return filename

def load_density_map(lambda_z_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    """Загрузка карты плотности из файла"""
    filename = get_cache_filename(lambda_z_range, incl, plane_rotation, cache_dir)
    if os.path.exists(filename):
        print(f"  Загрузка из кэша: {filename}")
        return np.load(filename)
    return None

def orbit_map_original(orb_groups, incl):
    """Создание карты плотности через Agama-проекцию.

    Применяет полную Эйлерову матрицу agama.makeRotationMatrix(alpha, beta, gamma)
    — точно ту же, что Agama использует при построении LOSVD внутри Target.
    Это гарантирует, что плотность находится в той же плоскости неба,
    что и кинематика модели и наблюдения.

    Дополнительный параметр `incl` поворачивает камеру вокруг линии узлов
    (X-оси наблюдателя) — позволяет показать face-on (incl=0) или edge-on
    (incl=90) view, не пересчитывая исходную модель.
    """
    image = np.zeros((scale, scale))

    # Базовая матрица проекции из galactic frame в плоскость неба.
    # Берём первые две строки 3×3 матрицы — это X и Y в observer frame.
    M_agama = np.array(agama.makeRotationMatrix(
        alpha=np.radians(agama_alpha),
        beta=np.radians(agama_beta),
        gamma=np.radians(agama_gamma)))  # 3×3, galactic → observer

    # Дополнительный поворот камеры вокруг X-оси наблюдателя (для col 1
    # — чтобы перейти от текущего displayed_inclination к incl).
    # incl_extra = incl - displayed_inclination_col0
    # Но проще: применяем rotation_matrix_x(incl - real_inclination)
    # как доворот после Agama-проекции.
    # ВНИМАНИЕ: agama.makeRotationMatrix(beta) крутит вокруг X в ПРОТИВОПОЛОЖНУЮ
    # сторону к rotation_matrix_x → знак incl_extra инвертирован, иначе col1
    # (incl=90) даёт почти face-on вместо edge-on. col0 (incl_extra=0) не затронут.
    incl_extra = real_inclination - incl
    M_extra = rotation_matrix_x(np.radians(incl_extra))  # 3×3

    # Полная матрица проекции: сначала доворот в galactic frame
    # (rotation_matrix_x на incl_extra), потом Agama-проекция на небо.
    # Доворот применяется ДО Agama, потому что в galactic frame
    # большая ось диска совпадает с X — поворот вокруг X есть
    # поворот вокруг большой оси → правильный edge-on view.
    # M[0:2] @ position даст 2D-координаты в плоскости неба.
    M = (M_agama @ M_extra)[0:2]  # 2×3

    # Азимутальное усреднение (модель осесимметрична): размазка траектории по φ,
    # иначе редконаселённый-по-массе противовращ. диск даёт «шестерёнку».
    # Векторизовано по точкам.
    N_AZ = 24
    dphi = np.linspace(0.0, 2.0*np.pi, N_AZ, endpoint=False)
    cos_d = np.cos(dphi)[None, :]
    sin_d = np.sin(dphi)[None, :]
    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[model_index][orbi]
            try:
                with open(f"C:/Users/Кель Ыр/Desktop/Галактики/Статья-2/LEDA/results_C/orbits_M1e+07_O0_Rh30.6_Vh125_i27_a0_N2000/orbit_{orbi}.txt", "r") as f:
                    pts = np.array([t[:3] for t in (ln.split() for ln in f) if len(t) >= 3],
                                   dtype=float)
            except (FileNotFoundError, OSError, ValueError):
                continue
            if pts.ndim != 2 or pts.shape[0] == 0:
                continue
            R   = np.hypot(pts[:, 0], pts[:, 1])[:, None]
            ph  = np.arctan2(pts[:, 1], pts[:, 0])[:, None]
            cph = np.cos(ph); sph = np.sin(ph)
            xx = (R * (cph*cos_d - sph*sin_d)).ravel()
            yy = (R * (sph*cos_d + cph*sin_d)).ravel()
            zz = np.repeat(pts[:, 2], N_AZ)
            sky = M @ np.vstack([xx, yy, zz])
            xi = (sky[0] + scale/2).astype(int)
            yi = (sky[1] + scale/2).astype(int)
            ok = (xi >= 0) & (xi < scale) & (yi >= 0) & (yi < scale)
            np.add.at(image, (yi[ok], xi[ok]), weight/1000.0/N_AZ)

    image[image == 0] = np.nan
    return image

def orbit_map_rotated_plane(orb_groups, incl, plane_rotation):
    """Создание карты плотности с вращением плоскости проекции - черный фон для NaN"""
    image = np.zeros((scale, scale))
    
    rot_matrix = rotation_matrix_2d(-plane_rotation)  # знак инвертирован для согласования с кинематикой
    
    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[model_index][orbi]
            try:
                # ИЗМЕНЕНО: Папка с орбитами из второго кода
                with open(f"C:/Users/Кель Ыр/Desktop/Галактики/Статья-2/LEDA/results_C/orbits_M1e+07_O0_Rh30.6_Vh125_i27_a0_N2000/orbit_{orbi}.txt", "r") as f:
                    for str_point in f:
                        point = [float(i) for i in str_point.split(" ")]
                        if len(point) < 3:
                            continue
                            
                        position = np.array([point[0], point[1], point[2]])
                        position = position @ rotation_matrix_x(np.radians(incl))
                        
                        # ИЗМЕНЕНО: проекция на плоскость неба (X, Y) вместо
                        # исторической (X, Z). См. orbit_map_original выше.
                        proj_coords = np.array([position[0], position[1]])
                        rotated_coords = rot_matrix @ proj_coords
                        
                        x_ = int(rotated_coords[0] + scale/2)
                        y_ = int(rotated_coords[1] + scale/2)
                        
                        if 0 <= x_ < scale and 0 <= y_ < scale:
                            image[y_][x_] += weight/1000
            except FileNotFoundError:
                continue
    
    image[image == 0] = np.nan
    return image

def create_and_cache_density_map(lambda_z_range, incl, plane_rotation=0, use_rotated_plane=False):
    """
    Создание или загрузка карты плотности с кэшированием
    lambda_z_range: диапазон индексов λ_z (например, (13, 21) для со-вращения)
    """
    cache_dir = CACHE_DIR_ROTATED if use_rotated_plane else CACHE_DIR_ORIGINAL
    
    cached = load_density_map(lambda_z_range, incl, plane_rotation, cache_dir)
    if cached is not None:
        return cached
    
    print(f"  Создание новой карты плотности...")
    print(f"    λ_z индексы: {lambda_z_range}, incl: {incl}, plane_rotation: {plane_rotation}")
    print(f"    Метод: {'вращение плоскости проекции' if use_rotated_plane else 'оригинальный'}")
    
    # ИЗМЕНЕНО: Берем орбиты из массива 21x21 (радиус x λ_z)
    orb_list = [sorted_orbs[radius][lz] for lz in range(*lambda_z_range) for radius in range(21)]
    
    if use_rotated_plane:
        image = orbit_map_rotated_plane(orb_list, incl, plane_rotation)
    else:
        image = orbit_map_original(orb_list, incl)
    
    save_density_map(image, lambda_z_range, incl, plane_rotation, cache_dir)
    
    return image

# ========== КЛАСС ДЛЯ СОЗДАНИЯ ГРАФИКОВ ==========
class FormattedPlotCreator:
    """Создает графики с правильным оформлением"""
    
    def __init__(self):
        self.gamma = gamma
        self.tr = transforms.Affine2D().rotate_deg_around(scale/2, scale/2, -self.gamma)

        # Поворот плотности модели на PA (gamma1 в углах Эйлера Agama).
        # Для LEDA 2220522: alpha=0, beta=27°, gamma1=40° → density_PA = 40.
        self.density_PA = float(density_PA)
        # Транспозиция плотности — та же операция, что применяется к
        # модельной кинематике через flip_diagonal=True.
        self.density_transpose = bool(density_transpose)
        print("[DIAG] self.gamma=%s°, transform=rotate(%s°, CCW); density_PA=%s°, density_transpose=%s, flip_x=%s, flip_y=%s" %
              (self.gamma, -self.gamma, self.density_PA, self.density_transpose, density_flip_x, density_flip_y))
        # Доп. отражения (по умолчанию False).
        self.density_flip_x = bool(density_flip_x)
        self.density_flip_y = bool(density_flip_y)

        flags = []
        if self.density_transpose: flags.append("transpose")
        if self.density_flip_x:    flags.append("flip_x")
        if self.density_flip_y:    flags.append("flip_y")
        if self.density_PA != 0:   flags.append(f"PA={self.density_PA}°")
        if flags:
            print(f"  Преобразования плотности модели: {', '.join(flags)}")

        self.norm_density = colors.LogNorm(vmin=Ph_min, vmax=Ph_max)
        self.norm_photometry = colors.LogNorm(vmin=Phot_min, vmax=np.max(phot['CCD3'].data))

        # Границы компонент по индексам бинов λ_z (схема А, 21 бин 0..20).
        # Берутся из глобальных констант counter_rot_end_idx и co_rot_start_idx
        # — см. начало файла. Чтобы попробовать другое разбиение, поменяйте
        # эти два числа в шапке и перезапустите скрипт.
        self.counter_rot_start_idx = 0
        self.counter_rot_end_idx   = counter_rot_end_idx
        self.spherical_start_idx   = counter_rot_end_idx
        self.spherical_end_idx     = co_rot_start_idx
        self.co_rot_start_idx      = co_rot_start_idx
        self.co_rot_end_idx        = 21

        # Сводная информация в консоль для отслеживания, какое разбиение
        # сейчас используется (особенно важно при переборе вариантов).
        cr_lo, cr_hi = bin_index_to_lambda_z_range(self.counter_rot_start_idx,
                                                    self.counter_rot_end_idx)
        sp_lo, sp_hi = bin_index_to_lambda_z_range(self.spherical_start_idx,
                                                    self.spherical_end_idx)
        co_lo, co_hi = bin_index_to_lambda_z_range(self.co_rot_start_idx,
                                                    self.co_rot_end_idx)
        print(f"  Разбиение на компоненты:")
        print(f"    counter-rot: бины [{self.counter_rot_start_idx:2d}, "
              f"{self.counter_rot_end_idx:2d})  ↔  λ_z ∈ [{cr_lo:+.2f}, {cr_hi:+.2f}]")
        print(f"    spheroid:    бины [{self.spherical_start_idx:2d}, "
              f"{self.spherical_end_idx:2d})  ↔  λ_z ∈ [{sp_lo:+.2f}, {sp_hi:+.2f}]")
        print(f"    co-rot:      бины [{self.co_rot_start_idx:2d}, "
              f"{self.co_rot_end_idx:2d})  ↔  λ_z ∈ [{co_lo:+.2f}, {co_hi:+.2f}]")

    
    def create_formatted_plots(self, use_rotated_plane=False, rotate_kinematics=True, 
                              flip_diagonal_before_transform=False, filename_suffix=""):
        """
        Создание всех графиков
        """
        fig, axs = plt.subplots(4, 4, sharex=False, sharey=False)
        fig.subplots_adjust(hspace=0)
        fig.subplots_adjust(wspace=-0.58)
        
        fig.delaxes(axs[0, 0])
        
        # Общие подписи (сохраняем русские, как в первом коде)
        fig.text(0.19, 0.5, 'угл. сек.', va='center', rotation='vertical', fontsize=7)
        fig.text(0.83, 0.5, 'угл. сек.', va='center', rotation='vertical', fontsize=7)
        fig.text(0.5, 0.07, '    угл. сек.', ha='center', fontsize=7)
        
        plane_rotation = gamma if use_rotated_plane else 0
        
        # ===== РЯД 1: Наблюдательные данные =====
        print("Ряд 1: Наблюдательные данные...")
        
        # Фотометрия (0,1)
        self._plot_photometry(axs[0, 1], fig)
        
        # Наблюдательная скорость (0,2)
        velocity_data = spectr[6].data - np.nanmean(spectr[6].data)
        self._plot_observational_velocity(axs[0, 2], velocity_data, fig, 
                                         rotate_kinematics=rotate_kinematics)
        
        # Наблюдательная дисперсия (0,3)
        sigma_data = spectr[8].data
        self._plot_observational_sigma(axs[0, 3], sigma_data, fig, 
                                      rotate_kinematics=rotate_kinematics)
        
        # ===== РЯД 2: Со-вращающийся диск =====
        print("Ряд 2: Со-вращающийся диск...")
        self._plot_component_row(axs, 1, 
                               lambda_z_range=(self.co_rot_start_idx, self.co_rot_end_idx),
                               bounds=list(bin_index_to_lambda_z_range(
                                   self.co_rot_start_idx, self.co_rot_end_idx)),
                               component_name="Со-вращ. диск", 
                               use_rotated_plane=use_rotated_plane, 
                               plane_rotation=plane_rotation,
                               rotate_kinematics=rotate_kinematics,
                               flip_diagonal=flip_diagonal_before_transform,
                               fig=fig)
        
        # ===== РЯД 3: Контр-вращающийся диск =====
        print("Ряд 3: Контр-вращающийся диск...")
        self._plot_component_row(axs, 2, 
                               lambda_z_range=(self.counter_rot_start_idx, self.counter_rot_end_idx),
                               bounds=list(bin_index_to_lambda_z_range(
                                   self.counter_rot_start_idx, self.counter_rot_end_idx)),
                               component_name="Контр-вращ. диск",
                               use_rotated_plane=use_rotated_plane, 
                               plane_rotation=plane_rotation,
                               rotate_kinematics=rotate_kinematics,
                               flip_diagonal=flip_diagonal_before_transform,
                               fig=fig)
        
        # ===== РЯД 4: Сферическая компонента =====
        print("Ряд 4: Сферическая компонента...")
        self._plot_component_row(axs, 3, 
                               lambda_z_range=(self.spherical_start_idx, self.spherical_end_idx),
                               bounds=list(bin_index_to_lambda_z_range(
                                   self.spherical_start_idx, self.spherical_end_idx)),
                               component_name="Сферическая",
                               use_rotated_plane=use_rotated_plane, 
                               plane_rotation=plane_rotation,
                               rotate_kinematics=rotate_kinematics,
                               flip_diagonal=flip_diagonal_before_transform,
                               fig=fig)
        
        filename = f"losvd_and_surface_density_LEDA_2220522{filename_suffix}.pdf"
        print(f"Сохранение в {filename}...")
        fig.savefig(filename, format="pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)
        print("Готово!")
    
    def _plot_photometry(self, ax, fig):
        """Фотометрия - УБРАНЫ ПОДПИСИ К ОСЯМ"""
        # ИЗМЕНЕНО: Используем CCD3 для LEDA 2220522
        image = phot['CCD3'].data
        image[(image <= 0)] = Phot_min
        
        pc = ax.imshow(image, cmap=cmap_2, norm=self.norm_photometry, origin='lower')
        
        # Настройки осей для фотометрии - тики остаются, но без подписей
        ax.yaxis.set_major_locator(tcr.NullLocator())
        ax.xaxis.set_major_locator(tcr.NullLocator())
        ax.yaxis.set_major_formatter(tcr.NullFormatter())
        ax.xaxis.set_major_formatter(tcr.NullFormatter())
        # ИЗМЕНЕНО: Offsets из второго кода (-8.5 вместо 2.5)
        ax.yaxis.set_major_locator(tcr.IndexLocator(6, offset=-7.5))
        ax.xaxis.set_major_locator(tcr.IndexLocator(6, offset=-7.5))
        ax.yaxis.set_minor_formatter(tcr.NullFormatter())
        ax.xaxis.set_minor_formatter(tcr.NullFormatter())
        # ИЗМЕНЕНО: Offsets из второго кода (-2.5 вместо 0.5)
        ax.yaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))
        ax.xaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))
        # ИЗМЕНЕНО: уменьшен размер шрифта с 8 до 6 (в 1.5 раза)
        ax.tick_params(direction="in", which="major", width=0.3, length=2, labelsize=6)
        ax.tick_params(direction="in", which="minor", width=0.15, length=1)
        
        # ИЗМЕНЕНО: Позиция текста поправлена для scale=42
         # НОВАЯ ПОЗИЦИЯ: Поток сверху слева - выше (ближе к верху)
        ax.annotate(r"$Log_{10}(Поток)[отсчет/пиксель]$", (3, 92),  
                   fontsize=4, style="italic", path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        
        # Колорбар
        divider = make_axes_locatable(ax)
        ax_cb = divider.append_axes("top", size="5%", pad=-0.045)
        cb_ph = fig.colorbar(pc, ax=ax, cax=ax_cb, orientation='horizontal', location='top', 
                     format=tcr.NullFormatter(), ticks=[np.max(image)])
        ax.set_title(f'{np.log10(np.max(image)):.1f} ', loc='right', fontdict=fontdict, pad=0.5)
        ax.set_title(f'{np.log10(np.min(image)):.0f} ', loc='left', fontdict=fontdict, pad=0.5)
        # ИЗМЕНЕНО: Убраны вторичные оси с подписями для фотометрии
        # Больше не создаем secondary axes для фотометрии
    
    def _plot_observational_velocity(self, ax, data, fig, rotate_kinematics=True):
        """Наблюдательная скорость - ПОДПИСИ К ТИКАМ НА ВСЕХ Y-ОСЯХ ПОВЕРНУТЫ"""
        divider = make_axes_locatable(ax)
        ax_cb = divider.append_axes("top", size="5%", pad=-0.045)
        
        if rotate_kinematics:
            transform = self.tr + ax.transData
        else:
            transform = ax.transData
        
        pc = ax.imshow(data, cmap=cmap, transform=transform, 
                      vmin=Vel_min, vmax=Vel_max, origin='lower')
        
        ax.yaxis.set_major_locator(tcr.NullLocator())
        ax.xaxis.set_major_locator(tcr.NullLocator())
        fig.colorbar(pc, ax=ax, cax=ax_cb, orientation='horizontal', location='top', 
                    ticks=[Vel_min, Vel_max], format=tcr.NullFormatter())
        ax.set_title(f' {Vel_min}', loc='left', fontdict=fontdict, pad=0.5)
        ax.set_title(f'{Vel_max} ', loc='right', fontdict=fontdict, pad=0.5)
        
        # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,26)
        ax.annotate(r" $V_{0}[км/с]$", (1, 26), fontsize=6, style="italic")
        
        ax.yaxis.set_major_formatter(tcr.NullFormatter())
        ax.xaxis.set_major_formatter(tcr.NullFormatter())
        # ИЗМЕНЕНО: Offsets из второго кода (-3.5 вместо -7.5)
        ax.yaxis.set_major_locator(tcr.IndexLocator(6, offset=-7.5))
        ax.xaxis.set_major_locator(tcr.IndexLocator(6, offset=-7.5))
        ax.yaxis.set_minor_formatter(tcr.NullFormatter())
        ax.xaxis.set_minor_formatter(tcr.NullFormatter())
        # ИЗМЕНЕНО: Offsets из второго кода (0.5)
        ax.yaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))
        ax.xaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))
        # ИЗМЕНЕНО: уменьшен размер шрифта с 8 до 6
        ax.tick_params(direction="in", which="major", width=0.3, length=2, labelsize=6)
        ax.tick_params(direction="in", which="minor", width=0.15, length=1)
        
        # ИЗМЕНЕНО: Добавлены вторичные оси с арксеками
        secondary_Axis_X = ax.secondary_xaxis("top")
        secondary_Axis_Y = ax.secondary_yaxis("right")
        
        # Устанавливаем те же тики, что и на основных осях
        major_locs = np.arange(-6, 7, 3) * 2 + scale/2  # Из второго кода: 12 шаг
        secondary_Axis_X.set_xticks(major_locs)
        secondary_Axis_Y.set_yticks(major_locs)
        
        # ИЗМЕНЕНО: Убраны подписи к тикам на верхней оси (NullFormatter)
        # Подписи на правой оси остаются, но будут повернуты
        secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
        secondary_Axis_Y.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
        
        # ИЗМЕНЕНО: Тики в одинаковом формате как на основных осях
        # Правые оси (Y) - подписи повернуты на 90 градусов
        secondary_Axis_X.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6)
        secondary_Axis_Y.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6, rotation=90)
    
    def _plot_observational_sigma(self, ax, data, fig, rotate_kinematics=True):
        """Наблюдательная дисперсия - ПОДПИСИ К ТИКАМ НА ВСЕХ Y-ОСЯХ ПОВЕРНУТЫ"""
        divider = make_axes_locatable(ax)
        ax_cb = divider.append_axes("top", size="5%", pad=-0.045)
        ax_cb.yaxis.set_major_locator(tcr.NullLocator())
        
        if rotate_kinematics:
            transform = self.tr + ax.transData
        else:
            transform = ax.transData
        
        pc = ax.imshow(data, cmap=cmap, transform=transform, 
                      vmin=Sigma_min, vmax=Sigma_max, origin='lower')
        
        fig.colorbar(pc, ax=ax, cax=ax_cb, orientation='horizontal', location='top', 
                    ticks=[Sigma_min, Sigma_max], format=tcr.NullFormatter())
        ax.set_title(f' {Sigma_min}', loc='left', fontdict=fontdict, pad=0.5)
        ax.set_title(f'{Sigma_max} ', loc='right', fontdict=fontdict, pad=0.5)
        
        # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,26)
        ax.annotate(r" $\sigma[км/с]$", (1, 26), fontsize=6, style="italic")
        
        ax.yaxis.set_major_formatter(tcr.NullFormatter())
        ax.xaxis.set_major_formatter(tcr.NullFormatter())
        # ИЗМЕНЕНО: Offsets из второго кода (-3.5 вместо -7.5)
        ax.yaxis.set_major_locator(tcr.IndexLocator(6, offset=-7.5))
        ax.xaxis.set_major_locator(tcr.IndexLocator(6, offset=-7.5))
        ax.yaxis.set_minor_formatter(tcr.NullFormatter())
        ax.xaxis.set_minor_formatter(tcr.NullFormatter())
        # ИЗМЕНЕНО: Offsets из второго кода (0.5)
        ax.yaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))
        ax.xaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))
        # ИЗМЕНЕНО: уменьшен размер шрифта с 8 до 6
        ax.tick_params(direction="in", which="major", width=0.3, length=2, labelsize=6)
        ax.tick_params(direction="in", which="minor", width=0.15, length=1)
        
        # ИЗМЕНЕНО: Добавлены вторичные оси с арксеками
        secondary_Axis_X = ax.secondary_xaxis("top")
        secondary_Axis_Y = ax.secondary_yaxis("right")
        
        # Устанавливаем те же тики, что и на основных осях
        major_locs = np.arange(-6, 7, 3) * 2 + scale/2  # Из второго кода: 12 шаг
        secondary_Axis_X.set_xticks(major_locs)
        secondary_Axis_Y.set_yticks(major_locs)
        
        # ИЗМЕНЕНО: Убраны подписи к тикам на верхней оси (NullFormatter)
        # Подписи на правой оси остаются, но будут повернуты
        secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
        secondary_Axis_Y.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
        
        # ИЗМЕНЕНО: Тики в одинаковом формате как на основных осях
        # Правые оси (Y) - подписи повернуты на 90 градусов
        secondary_Axis_X.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6)
        secondary_Axis_Y.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6, rotation=90)
    
    def _plot_component_row(self, axs, row, lambda_z_range, bounds, component_name, 
                           use_rotated_plane=False, plane_rotation=0, 
                           rotate_kinematics=True, flip_diagonal=False, fig=None):
        """Построение строки для одного компонента - ИСПРАВЛЕННЫЕ ПОЗИЦИИ НАДПИСЕЙ"""
        # Поворот плотности на density_PA реализуется через
        # orbit_map_rotated_plane (точный поворот точек траекторий ДО
        # дискретизации в пиксели). Это аналогично use_rotated_plane=True
        # в neworbmap.py — никакой scipy-интерполяции, никакой потери
        # масштаба и пикселей.
        # Если пользователь явно задал use_rotated_plane=True (например,
        # в creator2 с другим plane_rotation), его выбор имеет приоритет.
        effective_plane_rotation = plane_rotation
        effective_use_rotated_plane = use_rotated_plane
        if self.density_PA != 0 and not use_rotated_plane:
            effective_plane_rotation = self.density_PA
            effective_use_rotated_plane = True

        print(f"  Создание/загрузка плотности col 0 (i={incl_col0}°)...")
        lz_start_idx, lz_end_idx = lambda_z_range
        density_face = create_and_cache_density_map(
            lambda_z_range, incl_col0, effective_plane_rotation, effective_use_rotated_plane
        )
        
        print(f"  Создание/загрузка плотности col 1 (i={incl_col1}°)...")
        density_incl = create_and_cache_density_map(
            lambda_z_range, incl, effective_plane_rotation, effective_use_rotated_plane
        )

        # Преобразования плотности модели для согласования с системой
        # кинематики и наблюдений. См. константы density_* в начале файла.
        # Поворот уже применён выше через orbit_map_rotated_plane,
        # поэтому scipy.ndimage.rotate здесь не вызывается.
        if self.density_transpose:
            density_face = density_face.T
            density_incl = density_incl.T
        if self.density_flip_x:
            density_face = density_face[:, ::-1]
            density_incl = density_incl[:, ::-1]
        if self.density_flip_y:
            density_face = density_face[::-1, :]
            density_incl = density_incl[::-1, :]

        component_mass = 0.0
        for lz in range(lz_start_idx, lz_end_idx):
            for radius in range(21):  # 21 радиус
                for orbi in sorted_orbs[radius][lz]:
                    component_mass += weights[model_index][orbi]
        
        # Считаем общую массу ВСЕХ орбит (сумма весов по всем орбитам)
        total_mass = 0.0
        for lz in range(21):  # 21 λ_z
            for radius in range(21):  # 21 радиус
                for orbi in sorted_orbs[radius][lz]:
                    total_mass += weights[model_index][orbi]
        
        # ИЗМЕНЕНО: Массовая доля как отношение массы компонента к общей массе
        mass_fraction_face = (component_mass / total_mass) * 100
        
        # Для первой строки добавляем общий колорбар плотности
        if row == 1:
            ax_cb = plt.axes((0.2243, 0.689, 0.144, 0.01))
            
            # ИЗМЕНЕНО: Используем одну цветовую карту с черным фоном для обоих методов
            cmap_to_use = cmap_2 if not use_rotated_plane else cmap_density_rotated
            
            pc_face = axs[row, 0].imshow(density_face, cmap=cmap_to_use, norm=self.norm_density, origin='lower')
            cb_ph = fig.colorbar(pc_face, ax=axs[row, 0], cax=ax_cb, orientation='horizontal', 
                               location='top', format=tcr.NullFormatter(), ticks=[Ph_min,Ph_max])
            #cb_ph.ax.invert_xaxis()
            # ИЗМЕНЕНО: Подпись как во втором коде
            axs[row, 0].set_title(f" {np.log10(Ph_max * M_total * ML):.1f}" , 
                                loc='right', fontdict=fontdict, pad=5)
            axs[row, 0].set_title(f" {np.log10(Ph_min * M_total * ML):.1f}" , 
                                loc='left', fontdict=fontdict, pad=5)
            #axs[row, 0].set_title(f'{0}  ', loc='right', fontdict=fontdict, pad=5)
            axs[2, 0].set_title(f"LEDA 2220522", fontdict=fontdict, pad=100)  # Из второго кода
            
            # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,28)
            axs[row, 0].annotate(r"$Log_{10}(\mu_{b})[M_{\odot}/пиксель]$", (1, 28),
                               fontsize=5, style="italic", path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        else:
            cmap_to_use = cmap_2 if not use_rotated_plane else cmap_density_rotated
            pc_face = axs[row, 0].imshow(density_face, cmap=cmap_to_use, norm=self.norm_density, origin='lower')
        
        # Плотность face-on (слева)
        self._configure_density_axes(axs[row, 0], row, col=0, show_title=(row==1))
        
        # ИСПРАВЛЕННАЯ СТРУКТУРА НАДПИСЕЙ для первого столбца:
        if row == 1:
            # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,23)
            axs[row, 0].annotate(r"$i =$" + f"{displayed_inclination_col0}" + r"$\degree$", (1, 23), fontsize=6,
                               path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        else:
            # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,28)
            axs[row, 0].annotate(r"$i =$" + f"{displayed_inclination_col0}" + r"$\degree$", (1, 28), fontsize=6,
                               path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        
        # Отношение масс снизу слева - немного выше
        # ИЗМЕНЕНО: Используем M_total вместо M_Σ как во втором коде
        # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,1)
        axs[row, 0].annotate(r"$M_{c}/M_{total} =$" + f"{int(mass_fraction_face)}" + "%", 
                           (1, 1), fontsize=6, path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        
        # Плотность inclined (справа от face-on)
        cmap_to_use = cmap_2 if not use_rotated_plane else cmap_density_rotated
        pc_incl = axs[row, 1].imshow(density_incl, cmap=cmap_to_use, norm=self.norm_density, origin='lower')
        self._configure_density_axes(axs[row, 1], row, col=1)
        
        # ИСПРАВЛЕННАЯ СТРУКТУРА для второго столбца: i сверху слева
        # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,28)
        axs[row, 1].annotate(r"$i =$" + f"{displayed_inclination_col1}" + r"$\degree$", 
                           (1, 28), fontsize=6, style="italic",
                           path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        
        # λ_z снизу слева - немного выше
        # ИЗМЕНЕНО: Правильное оформление λ с сабскриптом z
        # ИЗМЕНЕНО: Позиция поправлена для scale=42 (было 1,1)
        lambda_z_text = f"{bounds[1]:.2f} > λ$_z$ > {bounds[0]:.2f}"
        axs[row, 1].annotate(lambda_z_text, (1, 1), fontsize=6,
                           path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        
        # Кинематические карты
        print(f"  Создание кинематических карт...")
        # ИЗМЕНЕНО: Используем функцию с индексами напрямую
        lz_start_idx, lz_end_idx = lambda_z_range
        kin_maps = DYN_COMP_LOSVD_MAP_INDICES(lz_start_idx, lz_end_idx)
        
        # Кинематика: скорость
        velocity_map = kin_maps[0]
        if flip_diagonal:
            velocity_map = velocity_map.T
        
        if rotate_kinematics:
            velocity_transform = self.tr + axs[row, 2].transData
        else:
            velocity_transform = axs[row, 2].transData
        
        velocity_map_for_display = velocity_map.copy()
        velocity_map_for_display[np.isnan(velocity_map_for_display)] = Vel_min - 1
        
        pc_velocity = axs[row, 2].imshow(velocity_map_for_display, cmap=cmap, 
                                        transform=velocity_transform,
                                        vmin=Vel_min, vmax=Vel_max,
                                        origin='lower')
        self._configure_kinematic_axes(axs[row, 2], row, col=2, is_sigma=False)
        
        # Кинематика: дисперсия
        sigma_map = kin_maps[1]
        if flip_diagonal:
            sigma_map = sigma_map.T
        
        if rotate_kinematics:
            sigma_transform = self.tr + axs[row, 3].transData
        else:
            sigma_transform = axs[row, 3].transData
        
        sigma_map_for_display = sigma_map.copy()
        sigma_map_for_display[np.isnan(sigma_map_for_display)] = Sigma_min - 1
        
        pc_sigma = axs[row, 3].imshow(sigma_map_for_display, cmap=cmap, 
                                     transform=sigma_transform,
                                     vmin=Sigma_min, vmax=Sigma_max,
                                     origin='lower')
        self._configure_kinematic_axes(axs[row, 3], row, col=3, is_sigma=True)
    
    def _configure_density_axes(self, ax, row, col, show_title=False):
        """Настройка осей для карт плотности - ПОДПИСИ К ТИКАМ НА ВСЕХ Y-ОСЯХ ПОВЕРНУТЫ"""
        ax.yaxis.set_major_locator(tcr.NullLocator())
        ax.xaxis.set_major_locator(tcr.NullLocator())
        
        # ИЗМЕНЕНО: Offsets из второго кода (-2.5 вместо -7.5)
        offset = -7.5
        
        if col == 0:
            if row == 3:
                ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
                ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
            else:
                ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
                ax.xaxis.set_major_formatter(tcr.NullFormatter())
        else:
            if row == 3:
                ax.yaxis.set_major_formatter(tcr.NullFormatter())
                ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
            else:
                ax.yaxis.set_major_formatter(tcr.NullFormatter())
                ax.xaxis.set_major_formatter(tcr.NullFormatter())
        
        ax.yaxis.set_major_locator(tcr.IndexLocator(6, offset=offset))
        ax.xaxis.set_major_locator(tcr.IndexLocator(6, offset=offset))
        ax.yaxis.set_minor_formatter(tcr.NullFormatter())
        ax.xaxis.set_minor_formatter(tcr.NullFormatter())
        ax.yaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))  # Из второго кода
        ax.xaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))  # Из второго кода
        # ИЗМЕНЕНО: уменьшен размер шрифта с 8 до 6
        ax.tick_params(direction="in", which="major", width=0.3, length=2, labelsize=6)
        ax.tick_params(direction="in", which="minor", width=0.15, length=1)
        
        # ИЗМЕНЕНО: Добавлены вторичные оси с арксеками
        secondary_Axis_X = ax.secondary_xaxis("top")
        secondary_Axis_Y = ax.secondary_yaxis("right")
        
        # Устанавливаем те же тики, что и на основных осях
        major_locs = np.arange(-6, 7, 3) * 2 + scale/2  # Из второго кода: 12 шаг
        secondary_Axis_Y.set_yticks(major_locs)
        secondary_Axis_X.set_xticks(major_locs)
        
        # Определяем, какие оси должны показывать подписи
        if col == 0 and row == 1:
            # Для первого столбца первой строки: верхняя ось без подписей, правая с подписями
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
        elif col == 1 and row == 1:
            # Для второго столбца первой строки: обе оси без подписей
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.NullFormatter())
        elif col == 3 and row == 1:
            # Для четвертого столбца первой строки: верхняя ось без подписей, правая с подписями
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
        elif col == 3 and row == 2:
            # Для четвертого столбца второй строки: верхняя ось без подписей, правая с подписями
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
        elif col == 3 and row == 3:
            # Для четвертого столбца третьей строки: верхняя ось без подписей, правая с подписями
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
        else:
            # Для всех остальных: обе оси без подписей
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.NullFormatter())
        
        # ИЗМЕНЕНО: Тики в одинаковом формате как на основных осях
        # Правые оси (Y) - подписи повернуты на 90 градусов
        secondary_Axis_X.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6)
        secondary_Axis_Y.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6, rotation=90)
    
    def _configure_kinematic_axes(self, ax, row, col, is_sigma=False):
        """Настройка осей для кинематических карт - ПОДПИСИ К ТИКАМ НА ВСЕХ Y-ОСЯХ ПОВЕРНУТЫ"""
        ax.yaxis.set_major_locator(tcr.NullLocator())
        ax.xaxis.set_major_locator(tcr.NullLocator())
        
        # ИЗМЕНЕНО: Offsets из второго кода (-2.5 вместо -7.5)
        offset = -7.5
        
        if row == 3:
            if col == 2:
                ax.yaxis.set_major_formatter(tcr.NullFormatter())
                ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
            elif col == 3:
                ax.yaxis.set_major_formatter(tcr.NullFormatter())
                ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
        else:
            ax.yaxis.set_major_formatter(tcr.NullFormatter())
            ax.xaxis.set_major_formatter(tcr.NullFormatter())
        
        ax.yaxis.set_major_locator(tcr.IndexLocator(6, offset=offset))
        ax.xaxis.set_major_locator(tcr.IndexLocator(6, offset=offset))
        ax.yaxis.set_minor_formatter(tcr.NullFormatter())
        ax.xaxis.set_minor_formatter(tcr.NullFormatter())
        ax.yaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))  # Из второго кода
        ax.xaxis.set_minor_locator(tcr.IndexLocator(2, offset=0.5))  # Из второго кода
        # ИЗМЕНЕНО: уменьшен размер шрифта с 8 до 6
        ax.tick_params(direction="in", which="major", width=0.3, length=2, labelsize=6)
        ax.tick_params(direction="in", which="minor", width=0.15, length=1)
        
        # ИЗМЕНЕНО: Добавлены вторичные оси с арксеками
        secondary_Axis_X = ax.secondary_xaxis("top")
        secondary_Axis_Y = ax.secondary_yaxis("right")
        
        # Устанавливаем те же тики, что и на основных осях
        major_locs = np.arange(-6, 7, 3) * 2 + scale/2  # Из второго кода: 12 шаг
        secondary_Axis_Y.set_yticks(major_locs)
        secondary_Axis_X.set_xticks(major_locs)
        
        if is_sigma:
            # Для дисперсии: верхняя ось без подписей, правая с подписями
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
        else:
            # Для скорости: обе оси без подписей
            secondary_Axis_X.xaxis.set_major_formatter(tcr.NullFormatter())
            secondary_Axis_Y.yaxis.set_major_formatter(tcr.NullFormatter())
        
        # ИЗМЕНЕНО: Тики в одинаковом формате как на основных осях
        # Правые оси (Y) - подписи повернуты на 90 градусов
        secondary_Axis_X.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6)
        secondary_Axis_Y.tick_params(axis='both', direction="in", which="both", 
                                   width=0.3, length=2, labelsize=6, rotation=90)

# ========== ОСНОВНОЙ КОД ==========
def main():
    """Создает два набора графиков"""
    print("=" * 60)
    print("СОЗДАНИЕ ГРАФИКОВ С ИСПРАВЛЕННЫМИ ПОЗИЦИЯМИ НАДПИСЕЙ")
    print("ДАННЫЕ ИЗ ВТОРОГО КОДА")
    print(f"Кэш оригинальных карт: {CACHE_DIR_ORIGINAL}")
    print(f"Кэш карт с вращением плоскости: {CACHE_DIR_ROTATED}")
    print("=" * 60)
    
    # 1. Набор 1: Как в первом коде (оригинальный)
    print("\nНАБОР 1: Как в первом коде (оригинальный)")
    print("- Кинематика: вращается на gamma через transform")
    print("- Модельная кинематика: отражается по диагонали")
    print("- Плотность: не вращается (оригинальный метод)")
    print("- ИСПРАВЛЕННЫЕ ПОЗИЦИИ НАДПИСЕЙ")
    print("- Убраны подписи к осям в фотометрии")
    print("- Подписи к тикам на ВСЕХ Y-осях повернуты на 90°")
    print("- Пределы скорости: ±175 км/с (увеличено с ±170)")
    print("- Пределы дисперсии: 80-250 км/с (из второго кода)")
    print("- Уменьшены цифры на тиках (в 1.5 раза)")
    print("- Черный фон для плохих значений на картах плотности")
    print("-" * 40)
    
    creator1 = FormattedPlotCreator()
    creator1.create_formatted_plots(
        use_rotated_plane=False,
        rotate_kinematics=False,
        # rotate_kinematics=False для LEDA: gamma=-40° включает
        # transform-поворот на +40°, который применился бы и к
        # наблюдательной, и к модельной кинематике. Но плотность
        # уже повёрнута на +40° через orbit_map_rotated_plane
        # (см. density_PA в шапке), а наблюдения и модель в
        # bin_scheme-coords уже согласованы между собой через .T —
        # поэтому дополнительного transform к кинематике не нужно.
        # Соответствует neworbmap.py creator2.
        flip_diagonal_before_transform=True,
        filename_suffix="_original"
    )
    
    # 2. Набор 2: С вращением плоскости проекции
    print("\nНАБОР 2: С вращением плоскости проекции")
    print("- Кинематика: НЕ вращается")
    print(f"- Плотность: вращение плоскости проекции на {-gamma}° (новый метод)")
    print("- Черный фон для повернутой плотности")
    print("- ИСПРАВЛЕННЫЕ ПОЗИЦИИ НАДПИСЕЙ")
    print("- Убраны подписи к осям в фотометрии")
    print("- Подписи к тикам на ВСЕХ Y-осях повернуты на 90°")
    print("- Пределы скорости: ±175 км/с (увеличено с ±170)")
    print("- Пределы дисперсии: 80-250 км/с (из второго кода)")
    print("- Уменьшены цифры на тиках (в 1.5 раза)")
    print("- Черный фон для плохих значений на картах плотности")
    print("-" * 40)
    
    creator2 = FormattedPlotCreator()
    creator2.create_formatted_plots(
        use_rotated_plane=True,
        rotate_kinematics=False,
        flip_diagonal_before_transform=True,
        filename_suffix="_plane_rotation"
    )
    
    print("\n" + "=" * 60)
    print("ВСЕ НАБОРЫ ГРАФИКОВ СОЗДАНЫ УСПЕШНО!")
    print("1. plots_formatted_21r_original.pdf - как в первом коде:")
    print("   - Кинематика: вращение через transform + отражение по диагонали")
    print("   - Плотность: оригинальный метод (черный фон для NaN)")
    print("   - Фотометрия: УБРАНЫ ПОДПИСИ К ОСЯМ")
    print("   - Все подписи к тикам на Y-осях (И ЛЕВЫХ И ПРАВЫХ) ПОВЕРНУТЫ на 90°")
    print(f"   - Скорость: ±{_Vel_lim} км/с (автолимит), дисперсия: {Sigma_min}–{Sigma_max} км/с (автолимит)")
    print("   - Уменьшены цифры на тиках (в 1.5 раза)")
    print("   - 21 бин по λ_z, границы: ±0.35 для дисков, [-0.35, 0.35] для сферической")
    print("   - 21 бин по радиусу")
    print("   - Черный фон для плохих значений на всех картах плотности")
    print("2. plots_formatted_21r_plane_rotation.pdf - с вращением плоскости проекции:")
    print("   - Кинематика: НЕ вращается + отражение по диагонали")
    print("   - Плотность: вращение плоскости проекции (черный фон для NaN)")
    print("   - Фотометрия: УБРАНЫ ПОДПИСИ К ОСЯМ")
    print("   - Все подписи к тикам на Y-осях (И ЛЕВЫХ И ПРАВЫХ) ПОВЕРНУТЫ на 90°")
    print(f"   - Скорость: ±{_Vel_lim} км/с (автолимит), дисперсия: {Sigma_min}–{Sigma_max} км/с (автолимит)")
    print("   - Уменьшены цифры на тиках (в 1.5 раза)")
    print("   - 21 бин по λ_z, границы: ±0.35 для дисков, [-0.35, 0.35] для сферической")
    print("   - 21 бин по радиусу")
    print("   - Черный фон для плохих значений на всех картах плотности")
    print("")
    print("Кэшированные карты сохранены в двух папках:")
    print(f"  - {CACHE_DIR_ORIGINAL}: карты без вращения плоскости (черный фон)")
    print(f"  - {CACHE_DIR_ROTATED}: карты с вращением плоскости проекции (черный фон)")
    print("=" * 60)

# ========== ЗАПУСК ==========
if __name__ == "__main__":
    main()
