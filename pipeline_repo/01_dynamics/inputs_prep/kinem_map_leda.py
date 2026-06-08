import numpy as np
import agama as _agama
import matplotlib.pyplot as plt
from matplotlib import ticker as tcr, transforms
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ===== ПАРАМЕТРЫ ДЛЯ LEDA 2220522 =====
M_total = 1.5e9
Max_rad = 10
Filename = "M1e+07_O0_Rh131_Vh324_i27_a0_N40000_R0.00_GH_DensityCylindricalLinear.npz"
R_scale = 0.5296
Im_scale = 32
distance = 108643
gamma = 0
model_index = 8
galaxy_name = "LEDA 2220522"

# Загрузка данных
archive = np.load(Filename, allow_pickle=True, encoding='latin1')
LOSVD = archive['LOSVD']
Lambda_z = archive['MOD_Lambda_z']
circ2 = archive["circ2"]
Vel = archive["Velocity"]
Upsilon = archive["Upsilon"]
weights = archive["weights"][model_index]
Rmean = archive["Rmean"] * R_scale / (distance * np.pi / 648000)
cosincl = archive["cosincl"]
ic = archive["ic"]
inttime = archive["inttime"]
DYN_COMP_LOSVD = archive["DYN_COMP_LOSVD"][model_index]
Units = _agama.getUnits()

spectr = fits.open("results_8254-1902_vorb020_md19_ad-1_nmom4.fits")
bin_scheme = np.loadtxt("bins_LEDA_2220522_Damirs.txt")

def create_kinematic_map(Lambda_z_bounds, GH_moment):
    LOSVD_plt = np.zeros((263, 47))
    
    for r in range(21):
        for l_z in range(int((1 + Lambda_z_bounds[0])*10.5), 
                         int((1 + Lambda_z_bounds[1])*10.5)):
            LOSVD_plt += DYN_COMP_LOSVD[l_z][r]
    
    GH_moments = _agama.ghMoments(
        matrix=LOSVD_plt * Upsilon[model_index]**-0.5,
        gridv=np.linspace(-250, 250, 46) * Upsilon[model_index]**0.5,
        degree=2, ghorder=6
    )[:, (1, 2, 6, 7, 8, 9)]
    
    kinem_map = np.full((Im_scale, Im_scale), np.nan, dtype=float)
    for bin_i in range(len(bin_scheme)):
        x, y, idx = bin_scheme[bin_i]
        kinem_map[int(x * 2 + Im_scale/2)][int(y * 2 + Im_scale/2)] = \
            GH_moments[int(idx)][GH_moment]
    
    return kinem_map

def setup_axis(ax, show_y=False, show_x=False, is_phot=False, offset_x=6.5, offset_y=6.5):
    def ticker_x_formatter(x, pos):
        return f'{abs((x - Im_scale/2)/2):.0f}'
    
    def ticker_y_formatter(y, pos):
        if is_phot:
            return f'{abs((y - Im_scale/2)+4):.0f}'
        return f'{abs((y - Im_scale/2)/2):.0f}'
    
    if show_y:
        ax.yaxis.set_major_formatter(tcr.FuncFormatter(ticker_y_formatter))
        ax.yaxis.set_major_locator(tcr.IndexLocator(5 if is_phot else 10, offset=offset_y))
    else:
        ax.yaxis.set_major_locator(tcr.NullLocator())
        ax.yaxis.set_major_formatter(tcr.NullFormatter())
    
    if show_x:
        ax.xaxis.set_major_formatter(tcr.FuncFormatter(ticker_x_formatter))
        ax.xaxis.set_major_locator(tcr.IndexLocator(10, offset=offset_x))
    else:
        ax.xaxis.set_major_locator(tcr.NullLocator())
        ax.xaxis.set_major_formatter(tcr.NullFormatter())
    
    minor_step = 1 if is_phot else 2
    ax.yaxis.set_minor_locator(tcr.IndexLocator(minor_step, offset=0.5))
    ax.xaxis.set_minor_locator(tcr.IndexLocator(minor_step, offset=0.5))
    ax.yaxis.set_minor_formatter(tcr.NullFormatter())
    ax.xaxis.set_minor_formatter(tcr.NullFormatter())
    
    ax.tick_params(direction="in", which="both", labelsize=9)

def add_colorbar(fig, ax, pc, vmin, vmax, label_format='{:.0f}'):
    divider = make_axes_locatable(ax)
    ax_cb = divider.append_axes("right", size="5%", pad=-0.06)
    cbar = fig.colorbar(pc, ax=ax, cax=ax_cb, orientation='vertical')
    
    cbar.set_ticks([])
    
    ax_cb.text(1.8, 0.05, label_format.format(vmin), 
               transform=ax_cb.transAxes, 
               rotation=90, 
               va='bottom', ha='center', fontsize=10)
    
    ax_cb.text(1.8, 0.95, label_format.format(vmax), 
               transform=ax_cb.transAxes, 
               rotation=90, 
               va='top', ha='center', fontsize=10)
    
    return cbar

tr = transforms.Affine2D().rotate_deg_around(Im_scale/2, Im_scale/2, gamma)

def plot_comparison_maps():
    # Загрузка данных
    phot_obs = np.loadtxt('converted_fits_LEDA_2220522_z.txt')
    phot_model = np.loadtxt('model_LEDA_2220522_z.txt')
    phot_res = np.abs(phot_obs - phot_model)
    
    # Исправление: использование маскированных массивов вместо None
    Vel_map_obs = spectr[6].data - np.nanmean(spectr[6].data)
    Vel_map_obs_mask = np.ma.array(Vel_map_obs, mask=(Vel_map_obs == 0))
    
    Sig_map_obs = spectr[8].data
    Sig_map_obs_mask = np.ma.array(Sig_map_obs, mask=(Sig_map_obs == 0))
    
    # Границы цветовых шкал
    Phot_min = np.min([phot_obs.min(), phot_model.min()])
    Phot_max = np.max([phot_obs.max(), phot_model.max()])
    Vel_min, Vel_max = -45, 45
    Sig_min, Sig_max = 0, np.nanmax(Sig_map_obs_mask)
    
    Phot_res_min = phot_res.min()
    Phot_res_max = phot_res.max()
    Vel_res_min = 0
    Sig_res_min = 0
    
    # Создаем копии цветовых карт с белым цветом для NaN
    cmap_phot = plt.cm.inferno_r.copy()
    cmap_phot.set_bad('white')
    
    cmap_kinem = plt.cm.jet.copy()
    cmap_kinem.set_bad('white')
    
    # ===== ПАРАМЕТРЫ КОМПОНОВКИ =====
    fig_width = 8.5
    fig_height = 8
    
    left_margin = 0.08
    right_margin = 0.95
    
    left_area_width = 0.56
    gap_width = 0.01
    right_area_width = 0.33
    
    wspace_left = -0.06
    
    total_available = right_margin - left_margin
    left_area_end = left_margin + left_area_width * total_available
    gap_start = left_area_end
    gap_end = gap_start + gap_width * total_available
    right_area_start = gap_end
    
    # Создаем фигуру с белым фоном
    fig = plt.figure(figsize=(fig_width, fig_height), facecolor='white')
    
    from matplotlib.gridspec import GridSpec
    
    # Левая область: столбцы 1 и 2
    gs_left = GridSpec(3, 2, figure=fig,
                      left=left_margin,
                      right=left_area_end,
                      hspace=0,
                      wspace=wspace_left)
    
    # Правая область: столбец 3
    gs_right = GridSpec(3, 1, figure=fig,
                       left=right_area_start,
                       right=right_margin,
                       hspace=0)
    
    Lambda_z_bounds = [-1, 1]
    
    # === РЯД 1: Фотометрия ===
    ax1 = fig.add_subplot(gs_left[0, 0])
    pc1 = ax1.imshow(phot_obs, cmap=cmap_phot, transform=tr + ax1.transData,
                    origin='lower', vmin=Phot_min, vmax=Phot_max)
    setup_axis(ax1, show_y=True, is_phot=True, offset_y=2.5)
    ax1.annotate(r"$Log_{10}(Flux)[counts/pix.]$", (1, 23), 
                fontsize=7, style="italic")
    ax1.set_title("obs.", fontsize=19)
    
    ax2 = fig.add_subplot(gs_left[0, 1])
    pc2 = ax2.imshow(phot_model, cmap=cmap_phot, transform=tr + ax2.transData,
                    origin='lower', vmin=Phot_min, vmax=Phot_max)
    setup_axis(ax2, is_phot=True, offset_y=2.5)
    add_colorbar(fig, ax2, pc2, round(Phot_min), round(Phot_max))
    ax2.set_title("model", fontsize=19)
    
    ax3 = fig.add_subplot(gs_right[0, 0])
    pc3 = ax3.imshow(phot_res, cmap=cmap_phot, transform=tr + ax3.transData,
                    origin='lower', vmin=Phot_res_min, vmax=Phot_res_max)
    setup_axis(ax3, is_phot=True, offset_y=2.5)
    add_colorbar(fig, ax3, pc3, round(Phot_res_min), round(Phot_res_max))
    ax3.annotate(galaxy_name, (9, 23), fontsize=11, style="italic")  # Изменена позиция
    ax3.set_title("residual", fontsize=19)
    
    # === РЯД 2: Скорость ===
    ax4 = fig.add_subplot(gs_left[1, 0])
    pc4 = ax4.imshow(Vel_map_obs_mask, cmap=cmap_kinem, transform=tr + ax4.transData,
                    origin='lower', vmin=Vel_min, vmax=Vel_max)
    setup_axis(ax4, show_y=True)
    ax4.set_ylabel("arcsec", fontsize=11)
    ax4.annotate(r"$V_{0}[\rm{km/s}]$", (2, 28), fontsize=11, style="italic")
    
    ax5 = fig.add_subplot(gs_left[1, 1])
    vel_model = create_kinematic_map(Lambda_z_bounds, 0)
    pc5 = ax5.imshow(vel_model, cmap=cmap_kinem, transform=tr + ax5.transData,
                    origin='lower', vmin=Vel_min, vmax=Vel_max)
    setup_axis(ax5)
    add_colorbar(fig, ax5, pc5, Vel_min, Vel_max, label_format='{:.0f}')
    
    ax6 = fig.add_subplot(gs_right[1, 0])
    vel_res = np.abs(create_kinematic_map(Lambda_z_bounds, 0) - Vel_map_obs_mask)
    Vel_res_max = np.nanmax(vel_res)
    pc6 = ax6.imshow(vel_res, cmap=cmap_kinem, transform=tr + ax6.transData,
                    origin='lower', vmin=Vel_res_min, vmax=Vel_res_max)
    setup_axis(ax6)
    add_colorbar(fig, ax6, pc6, Vel_res_min, Vel_res_max, label_format='{:.0f}')
    
    # === РЯД 3: Дисперсия ===
    ax7 = fig.add_subplot(gs_left[2, 0])
    pc7 = ax7.imshow(Sig_map_obs_mask, cmap=cmap_kinem, transform=tr + ax7.transData,
                    origin='lower', vmin=Sig_min, vmax=Sig_max)
    setup_axis(ax7, show_y=True, show_x=True)
    ax7.annotate(r"$\sigma[\rm{km/s}]$", (2, 28), fontsize=11, style="italic")
    
    ax8 = fig.add_subplot(gs_left[2, 1])
    sig_model = create_kinematic_map(Lambda_z_bounds, 1)
    pc8 = ax8.imshow(sig_model, cmap=cmap_kinem, transform=tr + ax8.transData,
                    origin='lower', vmin=Sig_min, vmax=Sig_max)
    setup_axis(ax8, show_x=True)
    ax8.set_xlabel("arcsec", fontsize=11)
    add_colorbar(fig, ax8, pc8, round(Sig_min), round(Sig_max))
    
    ax9 = fig.add_subplot(gs_right[2, 0])
    sig_res = np.abs(create_kinematic_map(Lambda_z_bounds, 1) - Sig_map_obs_mask)
    Sig_res_max = np.nanmax(sig_res)
    pc9 = ax9.imshow(sig_res, cmap=cmap_kinem, transform=tr + ax9.transData,
                    origin='lower', vmin=Sig_res_min, vmax=Sig_res_max)
    setup_axis(ax9, show_x=True)
    add_colorbar(fig, ax9, pc9, round(Sig_res_min), round(Sig_res_max))
    
    # Сохранение в EPS с белым фоном
    fig.savefig("leda_2220522_model_vs_data.png", 
                format="png", 
                dpi=400, 
                bbox_inches='tight',
                pad_inches=0.05,
                facecolor='white',
                edgecolor='white',
                transparent=False)
    
    plt.close(fig)

# === ЗАПУСК ===
plot_comparison_maps()
