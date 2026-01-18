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


# ========== КОНСТАНТЫ И НАСТРОЙКИ ==========
scale = 42  # Из второго кода
arcsec_to_kpc = 0.5696  # Из второго кода
N = 40000
M_total = 4.5e9  # Из второго кода
distance = 117490000  # Из второго кода
incl = 48  # Из второго кода
model_index = 3  # Из второго кода
gamma = -270  # Химия отзеркалена по оси у

# Колормап для химии
cmap_chem = plt.cm.jet
cmap_chem.set_bad(color='white')

# Путь для кэширования карт плотности
CACHE_DIR_ORIGINAL = "density_cache_original_21r"
CACHE_DIR_ROTATED = "density_cache_rotated_21r"
for cache_dir in [CACHE_DIR_ORIGINAL, CACHE_DIR_ROTATED]:
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

# ========== ЗАГРУЗКА ДАННЫХ ==========
print("Загрузка данных...")
# ИЗМЕНЕНО: Имя файла из второго кода
Filename = "M1e+07_O0_Rh107_Vh154_i42_a0_N40000_R0.00_GH_DensitySphHarm.npz"
archive = np.load(Filename, allow_pickle=True, encoding='latin1')

weights = archive["weights"]
ML = archive["Upsilon"][model_index]
lambda_z_list = archive["MOD_Lambda_z"]
Rmean_list = archive["Rmean"]
DYN_COMP_LOSVD = archive["DYN_COMP_LOSVD"][model_index]

# ИЗМЕНЕНО: Файлы спектра, фотометрии и bin_scheme из второго кода
#spectr = fits.open("results_8992-3704_vorb020_md19_ad-1_nmom4.fits")
#phot = fits.open("mosaic-00126063-PGC_35706-z-CCD4-image.fits")
bin_scheme = np.loadtxt("bins_PGC35706_Damirs.txt")
chem = fits.open("results_8992-3704_vorb020_md19_ad-1_nmom4.fits")

# ========== ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ==========
def DYN_COMP_LOSVD_MAP_INDICES(lz_start_idx, lz_end_idx):
    """
    Создание карт кинематики для заданного диапазона индексов λ_z
    lz_start_idx: начальный индекс (включительно)
    lz_end_idx: конечный индекс (не включая)
    """
    LOSVD_plt = np.zeros((263,47))
    
    # ИЗМЕНЕНО: Суммируем по заданным индексам
    for r in range(21):
        for l_z in range(lz_start_idx, lz_end_idx):
            LOSVD_plt += DYN_COMP_LOSVD[l_z][r]
    
    # Перевод из splines в гермитовы моменты
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
            x = int(bin_scheme[bin_i][0] * 2 + scale/2)
            y = int(bin_scheme[bin_i][1] * 2 + scale/2)
            if 0 <= x < scale and 0 <= y < scale:
                KINEM_MAP[y][x] = GH_moments[int(bin_scheme[bin_i][2])][gh_moment]
        LOSVD_MAP.append(KINEM_MAP)
    
    return LOSVD_MAP

def lambda_z_to_index(lambda_z):
    """Преобразование λ_z в индекс (21 бин)"""
    return int((1 + lambda_z) * 10)

# ========== ОТРИСОВКА ОРБИТ С КЭШИРОВАНИЕМ ==========
print("Сортировка орбит...")
# ИЗМЕНЕНО: 21x21 (радиус x λ_z) - возвращаем оригинальную размерность
sorted_orbs = [[[] for _ in range(21)] for _ in range(21)]

dyn_comps_data_map = np.full((21,21), -1) # карта radius x lambda_z с индексами динамических компонентов
dyn_comps_data = []
R_max = np.mean(Rmean_list) * 3

for orbi in range(N):
    Rmean = Rmean_list[orbi]
    lambda_z = lambda_z_list[orbi]
    # ИЗМЕНЕНО: Умножаем на 20 (21 бин по радиусу)
    idx_r = int((Rmean/R_max)*20)
    idx_lz = int((1 + lambda_z) * 10)  # Индекс от 0 до 20 (21 бин по λ_z)
    
    if  0 <= idx_r <= 20 and 0 <= idx_lz <= 20 and len(sorted_orbs[idx_r][idx_lz]) == 0:
        dyn_comps_data_map[idx_r][idx_lz] = np.max(dyn_comps_data_map) + 1
        dyn_comps_data.append([idx_r,idx_lz])
        sorted_orbs[idx_r][idx_lz].append(orbi)
    
    elif 0 <= idx_r <= 20 and 0 <= idx_lz <= 20:
        sorted_orbs[idx_r][idx_lz].append(orbi)

dyn_comps_data = np.array(dyn_comps_data) # информация о radius и lambda_z каждого компонента

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

def get_cache_filename(lambda_z_range, radius_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    """Генерация имени файла для кэша"""
    lz_min, lz_max = lambda_z_range
    r_min, r_max = radius_range
    if plane_rotation == 0:
        return os.path.join(cache_dir, f"density_lz_{lz_min}_{lz_max}_r_{r_min}_{r_max}_incl_{incl}.npy")
    else:
        return os.path.join(cache_dir, f"density_lz_{lz_min}_{lz_max}_r_{r_min}_{r_max}_incl_{incl}_plane_rot_{plane_rotation}.npy")

def save_density_map(density_map, lambda_z_range, radius_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    """Сохранение карта плотности в файл"""
    filename = get_cache_filename(lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    np.save(filename, density_map)
    #print(f"  Карта сохранена в {filename}")
    return filename

def load_density_map(lambda_z_range, radius_range, incl, plane_rotation=0, cache_dir=CACHE_DIR_ORIGINAL):
    """Загрузка карты плотности из файла"""
    filename = get_cache_filename(lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    if os.path.exists(filename):
        #print(f"  Загрузка из кэша: {filename}")
        return np.load(filename)
    return None

def orbit_map_original(orb_groups, incl):
    """Оригинальная функция создания карты плотности - БЕЗ вращения плоскости проекции"""
    
    image = np.zeros((scale, scale)) 
    
    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[model_index][orbi]
            try:
                # ИЗМЕНЕНО: Папка с орбитами из второго кода
                with open(f"orbits_M1e+07_O0_Rh107_Vh154_i42_a0_N40000/orbit_{orbi}.txt", "r") as f:
                    for str_point in f:
                        point = [float(i) for i in str_point.split(" ")]
                        if len(point) < 3:
                            continue
                        
                        position = np.array([point[0], point[1], point[2]])
                        position = position @ rotation_matrix_x(np.radians(incl))
                        
                        x_ = int(position[0] + scale/2)
                        y_ = int(position[2] + scale/2)
                        
                        if 0 <= x_ < scale and 0 <= y_ < scale:
                            image[y_][x_] += weight/1000
            except FileNotFoundError:
                continue
    
    #image[image == 0] = np.nan
    return image

def orbit_map_rotated_plane(orb_groups, incl, plane_rotation):
    """Создание карты плотности с вращением плоскости проекции - черный фон для NaN"""
    
    image = np.zeros((scale, scale))
    
    rot_matrix = rotation_matrix_2d(plane_rotation)
    
    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[model_index][orbi]

            try:
                # ИЗМЕНЕНО: Папка с орбитами из второго кода
                with open(f"orbits_M1e+07_O0_Rh107_Vh154_i42_a0_N40000/orbit_{orbi}.txt", "r") as f:
                    for str_point in f:
                        point = [float(i) for i in str_point.split(" ")]
                        if len(point) < 3:
                            continue
                            
                        position = np.array([point[0], point[1], point[2]])
                        position = position @ rotation_matrix_x(np.radians(incl))
                        
                        proj_coords = np.array([position[0], position[2]])
                        rotated_coords = rot_matrix @ proj_coords
                        
                        x_ = int(rotated_coords[0] + scale/2)
                        y_ = int(rotated_coords[1] + scale/2)
                        
                        if 0 <= x_ < scale and 0 <= y_ < scale:
                            image[y_][x_] += weight/1000
            except FileNotFoundError:
                continue
    
    #image[image == 0] = np.nan
    return image

def create_and_cache_density_map(lambda_z_range, radius_range, incl, plane_rotation=0, use_rotated_plane=False):
    """
    Создание или загрузка карты плотности с кэшированием
    lambda_z_range: диапазон индексов λ_z (например, (13, 21) для со-вращения)
    """
    cache_dir = CACHE_DIR_ROTATED if use_rotated_plane else CACHE_DIR_ORIGINAL
    
    cached = load_density_map(lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    if cached is not None:
        return cached
    
    #print(f"  Создание новой карты плотности...")
    #print(f"    λ_z индексы: {lambda_z_range}, rad: {radius_range}, incl: {incl}, plane_rotation: {plane_rotation}")
    #print(f"    Метод: {'вращение плоскости проекции' if use_rotated_plane else 'оригинальный'}")
    
    # ИЗМЕНЕНО: Берем орбиты из массива 21x21 (радиус x λ_z)
    orb_list = [sorted_orbs[radius][lz] for lz in range(*lambda_z_range) for radius in range(*radius_range)]
    if use_rotated_plane:
        image = orbit_map_rotated_plane(orb_list, incl, plane_rotation)
    else:
        image = orbit_map_original(orb_list, incl)
    
    save_density_map(image, lambda_z_range, radius_range, incl, plane_rotation, cache_dir)
    
    return image

# ========== ПОПУЛЯЦИОННО - ДИНАМИЧЕСКОЕ МОДЕЛИРОВАНИЕ ==========

# Подготовка к созданию матрицы весов и создание векторов с наблюдаемой химией.

bin_map = chem[5].data  #
age_map = chem[22].data #
age_err = chem[23].data # shape = scale x scale
met_map = chem[24].data #
met_err = chem[25].data #

num_bin = np.max(bin_map) + 1
dyn_comp_num = len(dyn_comps_data)

bin_met = np.zeros(num_bin * 2 - 3) # Вектор наблюдаемой металличности. * 2 - 3 - вектор регуляризации
bin_age = np.zeros(num_bin * 2 - 3) # Вектор наблюдаемых возрастов. * 2 - 3 - вектор регуляризации

for x_ in range(0,scale):
    for y_ in range(0,scale):
        if bin_map[y_][x_] != -1:
            bin_met[bin_map[y_][x_]] = met_map[y_][x_]
            bin_age[bin_map[y_][x_]] = age_map[y_][x_]

def reg_matrix_1(dyn_comp_i, bin_i, dyn_comps_data_map, dyn_comps_data, lambda_reg):
    
    # ==========
    # Функция для создания матрицы регуляризации. возвращает 4 если dyn_comp_i = bin_i,
    # -1 если динамический компонент соседствует с dyn_comp_i = bin_i на карте radius x lambda_z
    # ==========
    
    if dyn_comp_i == bin_i: return 4 * lambda_reg
    try:
        if dyn_comp_i == dyn_comps_data_map[dyn_comps_data[bin_i][0] + 1][dyn_comps_data[bin_i][1] + 0]: return -1 * lambda_reg
    except: True
    try:
        if dyn_comp_i == dyn_comps_data_map[dyn_comps_data[bin_i][0] + 0][dyn_comps_data[bin_i][1] + 1]: return -1 * lambda_reg
    except: True
    try:
        if dyn_comp_i == dyn_comps_data_map[dyn_comps_data[bin_i][0] - 1][dyn_comps_data[bin_i][1] + 0]: return -1 * lambda_reg
    except: True
    try:
        if dyn_comp_i == dyn_comps_data_map[dyn_comps_data[bin_i][0] + 0][dyn_comps_data[bin_i][1] - 1]: return -1 * lambda_reg
    except: True
    return 0

# ==========

def create_weight_matrix_of_dyn_comps(lambda_reg):

    # ==========
    # Создание матрицы N_bins х N_comps (матрица весов) содержащую информацию о вкладе каждого компонента
    # в вес данного бина (вес компонента нормализуется на массу бина). Также к матрице весов добавляется
    # матрица регуляризации (N_bins - 3) х N_comps для преодоления деградации решения уравнения.
    # После создания матрица передаётся вместе с вектором наблюдаемых возрастов или металличности 
    # в scipy.optimize.lsq_linear.
    #
    # lambda_reg - "сила" регуляризации. чем больше значение, тем более "гладкое"(smooth) решение
    # ==========

    print("=" * 30)
    print("Создание матрицы весов для популяционно - динамического моделирования:")
    
    # Матрица весов(для удобства создаётся как N_comps х (N_bins - 3), позже траспонируется).  * 2 - 3 нужно для последующего добавления матрицы регуляризации
    weight_matrix_dyn_comps = np.zeros((dyn_comp_num, num_bin * 2 - 3)) 
    
    weight_map_binned = np.zeros((num_bin, )) # Вектор с массой каждого бина
    
    dyn_comp_ind = 0
    
    for t, lambda_z in zip(tqdm(range(0,21)), range(0,21)):
        for radius in range(0,21):
        
            if len(sorted_orbs[radius][lambda_z]) == 0: # Я не понимаю почему оно так работает(пропуск компонента если в нём нет орбит)
                continue
            
            
            dyn_comp_bins_weights = np.zeros((num_bin * 2 - 3 , )) # * 2 - 3 нужно для последующего добавления матрицы регуляризации
            lambda_z_range = (lambda_z, lambda_z + 1)
            radius_range = (radius, radius + 1)
            
            # Карта весов данного компонента
            dyn_comp_weight_map = create_and_cache_density_map(lambda_z_range,
            radius_range, incl, plane_rotation=gamma, 
            use_rotated_plane=True)
            
            # Назначение ненормализованного вклада компонента в каждый бин и общей массы бина
            for x_ in range(0,scale):
                for y_ in range(0,scale):
                    if bin_map[y_][x_] != -1:
                    
                        dyn_comp_bins_weights[bin_map[y_][x_]] += dyn_comp_weight_map[y_][x_]
                        weight_map_binned[bin_map[y_][x_]] += dyn_comp_weight_map[y_][x_]
            
            weight_matrix_dyn_comps[dyn_comp_ind] = dyn_comp_bins_weights
            dyn_comp_ind += 1
    

    
    print("=" * 30)
    print("Сделано.")
    print("Масса в матрице весов от общей массы галакткик: ", np.sum(weight_map_binned))
    print("Количество динамических компонентов: ", dyn_comp_num)
    print("Количество бинов: ", num_bin)
    
    #Нормализация весов

    for bin_i in range(0,num_bin):
        for dyn_comp_i in range(0,dyn_comp_num):
            if weight_matrix_dyn_comps[dyn_comp_i][bin_i] != 0:
                weight_matrix_dyn_comps[dyn_comp_i][bin_i] = weight_matrix_dyn_comps[dyn_comp_i][bin_i]/np.sum(weight_map_binned[bin_i])
    
    weight_matrix_dyn_comps = np.transpose(weight_matrix_dyn_comps)
        
    #Добавление матрицы регуляризации к матрице весов:

    for dyn_comp_i in range(0,dyn_comp_num):
        for bin_i in range(0,num_bin - 3):
        
            weight_matrix_dyn_comps[bin_i + num_bin][dyn_comp_i] = reg_matrix_1(dyn_comp_i, bin_i, dyn_comps_data_map, dyn_comps_data, lambda_reg)
    
    return weight_matrix_dyn_comps, weight_map_binned
    
# ========== PLOTS ==========
    
def tickers_radius_formatter(x, pos):
    return f'{(x / 21) * R_max:.1f}'

def tickers_lambda_z_formatter(y, pos):
    return f'{(y - 10)/10:.1f}'
    
def make_lambda_z_vs_radius_vs_chem_plot(model_chem, dyn_comps_data, chem_min, chem_max, chem = "chem"):

    image = np.full((21,21), np.nan)
    
    for dyn_comp, dyn_comp_i in zip(dyn_comps_data, range(0, len(dyn_comps_data))):
    
        image[dyn_comp[1]][dyn_comp[0]] = model_chem[dyn_comp_i]

    fig, ax = plt.subplots()
    pc = ax.imshow(image,cmap = cmap_chem, vmin = chem_min, vmax = chem_max)
    fig.colorbar(pc,ax=ax,label= {chem})
    ax.set_ylabel("Lambda_z")
    ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_lambda_z_formatter))
    ax.set_xlabel("radius[arcsec]")
    ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_radius_formatter))
    fig.savefig("lambda_z_vs_radius_vs_" + chem + ".pdf",format = "pdf")


def make_model_chem_plot(model_chem, weight_matrix, chem_min, chem_max, chem = "chem"):
    
    image = np.full((scale,scale), np.nan)
    
    for x_ in range(0,scale):
        for y_ in range(0,scale):
            if bin_map[y_][x_] != -1:
                image[y_][x_] = np.sum(weight_matrix[bin_map[y_][x_]] * model_chem)
    
    fig, ax = plt.subplots()
    pc = ax.imshow(image,cmap = cmap_chem, vmin = chem_min, vmax = chem_max)
    fig.colorbar(pc,ax=ax,label= {chem})
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    fig.savefig("model_" + chem + ".pdf",format = "pdf")


def make_obs_chem_plot(chem_map, chem_min, chem_max, chem = "chem"):

    fig, ax = plt.subplots()
    pc = ax.imshow(chem_map,cmap = cmap_chem, vmin = chem_min, vmax = chem_max)
    fig.colorbar(pc,ax=ax,label= {chem})
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    fig.savefig("obs_" + chem + ".pdf",format = "pdf")

def make_model_dencity_plot(weight_map_binned):
    
    image = np.full((scale,scale), np.nan)
    
    for x_ in range(0,scale):
        for y_ in range(0,scale):
            if bin_map[y_][x_] != -1:
                image[y_][x_] = weight_map_binned[bin_map[y_][x_]]
    
    fig, ax = plt.subplots()
    pc = ax.imshow(image,cmap = cmap_chem, norm = "log")
    fig.colorbar(pc,ax=ax,label= "weight")
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    fig.savefig("model_dencity" + ".pdf",format = "pdf")

# ========== ОСНОВНОЙ КОД ==========
def main():
    
    weight_matrix, weight_map_binned = create_weight_matrix_of_dyn_comps(lambda_reg = 0.1)
    
    age_min, age_max = 0, 14
    met_min, met_max = -2, 1
    
    print("=" * 30)
    print("Решение системы уравнений с помощью метода наименьших квадратов.")
    print("Используется scipy.optimize.lsq_linear, method='bvls':")
    
    phi_age = scipy.optimize.lsq_linear(weight_matrix, 
    bin_age, bounds=(age_min, age_max), method='bvls', 
    tol=1e-30, lsq_solver=None, lsmr_tol=None, 
    max_iter=None, verbose=1, lsmr_maxiter=None)
    
    phi_met = scipy.optimize.lsq_linear(weight_matrix, 
    bin_met, bounds=(met_min, met_max), method='bvls', 
    tol=1e-30, lsq_solver=None, lsmr_tol=None, 
    max_iter=None, verbose=1, lsmr_maxiter=None)
    
    print("=" * 30)
    print("Результат scipy.optimize.lsq_linear для возрастов:")
    print(phi_age)
    
    print("=" * 30)
    print("Результат scipy.optimize.lsq_linear для металличности:")
    print(phi_met)
    
    # ========== создание графиков на основе полученных данных ==========
    
    model_age = phi_age.x # Вектор содержащий информацию о возрасте каждого динамического компонента
    
    model_met = phi_met.x # Вектор содержащий информацию о металличности каждого динамического компонента
    
    make_lambda_z_vs_radius_vs_chem_plot(model_age,dyn_comps_data, age_min, age_max, chem = "age")
    
    make_lambda_z_vs_radius_vs_chem_plot(model_met,dyn_comps_data, met_min, met_max, chem = "met")
    
    make_model_chem_plot(model_age, weight_matrix, age_min, age_max, chem = "age")
    
    make_model_chem_plot(model_met, weight_matrix, met_min, met_max, chem = "met")
    
    make_obs_chem_plot(age_map, age_min, age_max, chem = "age")
    
    make_obs_chem_plot(met_map, met_min, met_max, chem = "met")
    
    make_model_dencity_plot(weight_map_binned)
    
# ========== ЗАПУСК ==========
if __name__ == "__main__":
    main()
