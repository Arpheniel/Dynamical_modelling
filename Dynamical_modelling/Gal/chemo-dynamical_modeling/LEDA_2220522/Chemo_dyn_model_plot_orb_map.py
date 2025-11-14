import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import agama
from time import sleep
from tqdm import tqdm
from decimal import Decimal
from matplotlib import ticker as tcr
from astropy.io import fits
from matplotlib import transforms
import matplotlib.colors as colors
import matplotlib.patheffects as pe

cmap = plt.cm.nipy_spectral
cmap.set_bad(color='White') # Set 'bad' values (NaNs) to white
cmap_2 = plt.cm.gist_heat
cmap_2.set_bad(color='White') # Set 'bad' values (NaNs) to white
cmap.set_under(color = "white")

np.set_printoptions(threshold=sys.maxsize)
plt.get_cmap('inferno')
scale = 32
kpc_to_arcsec = 1/0.56
N = 5000  # number of orbits
M_total = 1.5e9  # in units of M_sun
ML = 6 #M/L
distance = 108643000# in [pc]
incl = 63
Filename = "M1e+07_O0_Rh100_Vh190_i27_a0_N5000_R1.00_GH_DensityCylindricalLinear.npz"

archive = np.load(Filename, allow_pickle=True, encoding='latin1')

weights = archive["weights"]# * 6.75e10
lambda_z_list = archive["MOD_Lambda_z"]
Rmean_list = archive["Rmean"]
DYN_COMP_LOSVD = archive["DYN_COMP_LOSVD"][0]

spectr = fits.open("manga-8254-1902-MAPS-VOR10-MILESHC-MASTARSSP.fits")

phot = fits.open("mosaic-00157902-LEDA2220522-z-CCD3-image.fits")

bin_scheme = np.loadtxt("bins_LEDA2220522.txt")


def DYN_COMP_LOSVD_MAP(bounds): # bounds = [min l_z, max L_z]

    LOSVD_plt = np.zeros((292,47))
    
    for r in range(21):
        for l_z in range(int((1 + bounds[0])*10.5),int((1 + bounds[1])*10.5)): # 0 = -1, 11 = 0 , 21 = 1
            LOSVD_plt +=  DYN_COMP_LOSVD[l_z][r]
    
    GH_moments = agama.ghMoments(matrix=LOSVD_plt * ML**-0.5,gridv=np.linspace(-250, 250, 46) * ML**0.5, degree=2, ghorder=6)[:,(1,2,6,7,8,9)]
    
    LOSVD_MAP = []
    
    for gh_moment in range(6):
        KINEM_MAP = np.full((scale,scale),None,dtype=float)
        for bin_i in range(len(bin_scheme)):
            KINEM_MAP[int(bin_scheme[bin_i][0]+scale/2)][int(bin_scheme[bin_i][1]+scale/2)] = GH_moments[int(bin_scheme[bin_i][2])][gh_moment]
        LOSVD_MAP.append(KINEM_MAP)
    return LOSVD_MAP



sorted_orbs = []

for i in range(21):
    sorted_orbs.append([])
    for j in range(21):
        sorted_orbs[i].append([])

for orbi in range(N):
    Rmean = Rmean_list[orbi]
    lambda_z = lambda_z_list[orbi]
    # i / a < r < (i + 1) / a
    # i / b + min(lambda_z_list) < lambda_z < (i + 1) / b + min(lambda_z_list)
    if int((Rmean/(np.mean(Rmean_list)*3))*20) <= 20 and int((1 + lambda_z)*10) <= 20:
        sorted_orbs[int((Rmean/(np.mean(Rmean_list)*3))*20)][int((1 + lambda_z)*10)].append(orbi)


#print(sorted_orbs[[1,2,3,4]][[1,2,3,4]])
#print(sum(ar_dict.values()) / len(ar_dict.values()))
#print(min(lambda_z_list), max(lambda_z_list), min(ar_dict.values()), max(ar_dict.values()))


def rotation_matrix_x(a):
    return np.array([[1, 0, 0],
                     [0, np.cos(a), -np.sin(a)],
                     [0, np.sin(a), np.cos(a)]])


def rotation_matrix_y(a):
    return np.array([[np.cos(a), 0, np.sin(a)],
                     [0, 1, 0],
                     [-np.sin(a), 0, np.cos(a)]])


def rotation_matrix_z(a):
    return np.array([[np.cos(a), -np.sin(a), 0],
                     [np.sin(a), np.cos(a), 0],
                     [0, 0, 1]])

def tickers_X_formatter(x, pos):
    return f'{x/2 - 10.5}'

def tickers_Y_formatter(y, pos):
    return f'{y/2 - 10.5}' #Y is inverted

def orbit_map(orb_groups,incl):
    n = 0
    #width, height = 800, 800
    image = np.zeros((scale,scale)) 
    pixel_size = 1
    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[5][orbi]
            for str_point in open(f"orbits/orbit_{orbi}.txt"):
                #if weight >= 0.01: continue

                point = [float(i) for i in str_point.split(" ")]
                
                position = [point[0], point[1], point[2]]
                position = position @ rotation_matrix_x(np.radians(incl))
                position = position @ rotation_matrix_y(np.radians(0))
                position = position @ rotation_matrix_z(np.radians(0))
                
                x_ = int(position[0] *pixel_size*kpc_to_arcsec + scale/2*pixel_size)
                y_ = int(position[2] *pixel_size*kpc_to_arcsec + scale/2*pixel_size)
                try:
                    image[y_//pixel_size][x_//pixel_size] += weight/1000#/(distance * np.pi / 648000)*(ML*M_total)
                except:
                    pass
    print(np.sum(image))#/M_total)
    print(np.max(image))
    #fig, ax = plt.subplots()
    #pc = ax.imshow(image/1e7,norm="linear",cmap="cmap_2")
    #ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
    #ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
    #ax.yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
    #ax.xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
    #ax.set_xlabel("arcsec")
    #ax.set_ylabel("arcsec")
    #fig.colorbar(pc,ax=ax,label=r'$Mass, 10^6$ $M_{\odot}$',format = tcr.LogFormatter())
    #fig.savefig("density_map")
    return image


#orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(0,21) for radius in range(0,21)],42)

fig, axs = plt.subplots(4, 4, sharex = False, sharey = False)

#fig.subplots_adjust(left=-0.1, bottom=0.05, right=0.8, top=1, wspace=0, hspace=0)
fontdict = {'fontsize': 6}

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=-0.58)

def tickers_X_formatter(x, pos):
    #print(x,pos)
    return f'{int((x - scale/2)/2)}'
def tickers_Y_formatter(y, pos):
    #print(y,pos)
    return f'{int((y - scale/2)/2)}' #Y is inverted

#axs.ylabel("arcsec")
#axs.xlabel("arcsec")

#(14,21)
#(0,9)
#range(9,14)
#max_ph = 643669665

Ph_max = 0.024

norm=colors.LogNorm(vmin=0.000000001, vmax=Ph_max)

gamma = -45

Vel_max = 70
Vel_min = -70

Sigma_max = 100
Sigma_min = 0

tr = transforms.Affine2D().rotate_deg_around(scale/2,scale/2,gamma)

######### data ###########

image = spectr['STELLAR_VEL'].data
image[(image == 0)] = None
image[0][0] = -1000

divider = make_axes_locatable(axs[0,2])
ax_cb = divider.append_axes("top", size="5%", pad=-0.045)

pc = axs[0, 2].imshow(image,cmap=cmap,transform = tr + axs[0,2].transData,vmin = Vel_min, vmax = Vel_max)
axs[0, 2].yaxis.set_major_locator(tcr.NullLocator())
axs[0, 2].xaxis.set_major_locator(tcr.NullLocator())
fig.colorbar(pc, ax=axs[0, 2],cax = ax_cb,orientation = 'horizontal',location = 'top',ticks = [Vel_min,Vel_max],format = tcr.NullFormatter())
axs[0, 2].set_title(f' {Vel_min}',loc = 'left',fontdict = fontdict,pad = 0.5)
axs[0, 2].set_title(f'{Vel_max} ',loc = 'right',fontdict = fontdict,pad = 0.5)
axs[0,2].annotate(r" $V_{0}[km/s]$",(0,5),fontsize=5,style = "italic")
#axs[0,2].annotate(f"{Lambda_z_bound[1]} > $\lambda_z$ > {Lambda_z_bound[0]}",(25,39), fontsize=15)

axs[0, 2].yaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 2].xaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 2].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[0, 2].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[0, 2].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[0, 2].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[0, 2].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[0, 2].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[0, 2].tick_params(direction = "in",which = "major",width = 0.3,length = 2)
axs[0, 2].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)


image = spectr['STELLAR_SIGMA'].data
image[(image == 0)] = None
image[0][0] = -1000

divider = make_axes_locatable(axs[0,3])
ax_cb = divider.append_axes("top", size="5%", pad=-0.045)
ax_cb.yaxis.set_major_locator(tcr.NullLocator())

pc = axs[0, 3].imshow(image,cmap=cmap,transform = tr + axs[0,3].transData,vmin = Sigma_min, vmax = Sigma_max)

fig.colorbar(pc, ax=axs[0, 3],cax = ax_cb,orientation = 'horizontal',location = 'top',ticks = [Sigma_min,Sigma_max],format = tcr.NullFormatter())
axs[0, 3].set_title(f' {Sigma_min}',loc = 'left',fontdict = fontdict,pad = 0.5)
axs[0, 3].set_title(f'{Sigma_max} ',loc = 'right',fontdict = fontdict,pad = 0.5)
axs[0,3].annotate(r" $\sigma[km/s]$",(0,5),fontsize=5,style = "italic")

axs[0, 3].yaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 3].xaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 3].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[0, 3].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[0, 3].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[0, 3].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[0, 3].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[0, 3].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[0, 3].tick_params(direction = "in",which = "major",width = 0.3,length = 2)
axs[0, 3].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

fig.delaxes(axs[0,0])

image = phot['CCD3'].data
image[(image <= 0)] = 0.01

pc = axs[0, 1].imshow(image,cmap=cmap_2,norm = 'log')
#fig.colorbar(pc, ax=axs[0, 1] ,label='cos(incl)',orientation = 'horizontal',location = 'top')
axs[0, 1].yaxis.set_major_locator(tcr.NullLocator())
axs[0, 1].xaxis.set_major_locator(tcr.NullLocator())

axs[0, 1].yaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 1].xaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 1].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[0, 1].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[0, 1].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[0, 1].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[0, 1].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[0, 1].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

axs[0, 1].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[0, 1].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)


axs[0, 1].annotate(r"$M_{total} =" + f"{int(M_total * ML/10**9)}" + r"\cdot 10^{9} M_{\odot}$",(5,98), fontsize=3,path_effects=[pe.withStroke(linewidth=1, foreground="white")])

axs[0, 1].annotate(r"$Log_{10}(Flux)[counts/pixel]$",(5,10),fontsize=3,style = "italic",path_effects=[pe.withStroke(linewidth=1, foreground="white")])

ax_cb = plt.axes((0.361, 0.689, 0.0072, 0.191))

cb_ph = fig.colorbar(pc, ax=axs[0, 1],cax = ax_cb,orientation = 'vertical',location = 'left',format = tcr.NullFormatter(),ticks = [np.max(image)])

axs[0, 1].set_ylabel(f'{int(np.max(image))} ',loc = 'top',fontdict = fontdict,labelpad = 1)
######### co-rotating disk ###########

bounds = [0.3,1]

image = orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(14,21) for radius in range(0,21)],0)
#image[(image < 1)] = None

cmap_2.set_bad(color='White')

pc = axs[1, 0].imshow(image,cmap=cmap_2,norm=norm)
axs[1, 0].yaxis.set_major_locator(tcr.NullLocator())
axs[1, 0].xaxis.set_major_locator(tcr.NullLocator())

ax_cb = plt.axes((0.2243, 0.689, 0.144, 0.01))

cb_ph = fig.colorbar(pc, ax=axs[1, 0],cax = ax_cb,orientation = 'horizontal',location = 'top',format = tcr.NullFormatter(),ticks = [Ph_max])

cb_ph.ax.invert_xaxis()

axs[1, 0].set_title(f" {Ph_max * M_total * ML/10 ** 8}" + r"$\cdot 10^{8}$",loc = 'left',fontdict = fontdict,pad = 5)
axs[1, 0].set_title(f'{0}  ',loc = 'right',fontdict = fontdict,pad = 5)
axs[2, 0].set_title(f"LEDA 2220522",fontdict = fontdict,pad = 100)

axs[1, 0].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
axs[1, 0].xaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 0].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[1, 0].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[1, 0].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 0].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 0].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 0].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 0].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[1, 0].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

axs[1, 0].annotate(r"$Log_{10}(\mu_{b})[M_{\odot}/pixel]$",(1,3),fontsize=3,style = "italic",path_effects=[pe.withStroke(linewidth=1, foreground="white")])
axs[1, 0].annotate(r"$\varphi = 0\degree$",(25,3),fontsize=3,style = "italic",path_effects=[pe.withStroke(linewidth=1, foreground="white")])
axs[1, 0].annotate(r"$M_{c}/M_{total} =$" +  f"{int(np.sum(image) * 100)}" + "%",(1,30), fontsize=4,path_effects=[pe.withStroke(linewidth=1, foreground="white")])
axs[1, 0].annotate(f"{bounds[1]} > $\lambda_z$ > {bounds[0]-0.05}",(21,30), fontsize=3,path_effects=[pe.withStroke(linewidth=1, foreground="white")])

image = orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(14,21) for radius in range(0,21)],incl)
#image[(image < 1)] = None

pc = axs[1, 1].imshow(image,cmap=cmap_2,norm=norm)
axs[1, 1].yaxis.set_major_locator(tcr.NullLocator())
axs[1, 1].xaxis.set_major_locator(tcr.NullLocator())

axs[1, 1].yaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 1].xaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 1].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[1, 1].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[1, 1].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 1].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 1].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 1].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 1].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[1, 1].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

axs[1, 1].annotate(r"$\varphi =$" + f"{incl}" + r"$\degree$",(25,3),fontsize=3,style = "italic",path_effects=[pe.withStroke(linewidth=1, foreground="white")])

LOSVD_MAP = DYN_COMP_LOSVD_MAP(bounds)

image = LOSVD_MAP[0]
image[0][0] = -1000

pc = axs[1, 2].imshow(image,cmap=cmap,transform = tr + axs[1,2].transData,vmin = Vel_min, vmax = Vel_max)

axs[1, 2].yaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 2].xaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 2].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[1, 2].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[1, 2].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 2].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 2].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 2].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 2].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[1, 2].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

image = LOSVD_MAP[1]
image[0][0] = -1000

pc = axs[1, 3].imshow(image,cmap=cmap,transform = tr + axs[1,3].transData,vmin = Sigma_min, vmax = Sigma_max)
axs[1, 3].yaxis.set_major_locator(tcr.NullLocator())
axs[1, 3].xaxis.set_major_locator(tcr.NullLocator())

axs[1, 3].yaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 3].xaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 3].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[1, 3].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[1, 3].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 3].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[1, 3].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 3].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1, 3].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[1, 3].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

######### counter-rotating disk ###########

bounds = [-1,-0.3]

image = orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(0,9) for radius in range(0,21)],0)
#image[(image < 1)] = None

pc = axs[2, 0].imshow(image,cmap=cmap_2,norm=norm)
axs[2, 0].yaxis.set_major_locator(tcr.NullLocator())
axs[2, 0].xaxis.set_major_locator(tcr.NullLocator())

axs[2, 0].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
axs[2, 0].xaxis.set_major_formatter(tcr.NullFormatter())
axs[2, 0].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[2, 0].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[2, 0].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 0].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 0].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 0].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 0].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[2, 0].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

axs[2, 0].annotate(r"$M_{c}/M_{total} =$" +  f"{int(np.sum(image) * 100)}" + "%",(1,30), fontsize=4,path_effects=[pe.withStroke(linewidth=1, foreground="white")])

axs[2, 0].annotate(f"{bounds[1]+0.05} > $\lambda_z$ > {bounds[0]}",(20,30), fontsize=3,path_effects=[pe.withStroke(linewidth=1, foreground="white")])

image = orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(0,9) for radius in range(0,21)],incl)
#image[(image < 1)] = None

pc = axs[2, 1].imshow(image,cmap=cmap_2,norm=norm)
axs[2, 1].yaxis.set_major_locator(tcr.NullLocator())
axs[2, 1].xaxis.set_major_locator(tcr.NullLocator())

axs[2, 1].yaxis.set_major_formatter(tcr.NullFormatter())
axs[2, 1].xaxis.set_major_formatter(tcr.NullFormatter())
axs[2, 1].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[2, 1].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[2, 1].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 1].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 1].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 1].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 1].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[2, 1].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

LOSVD_MAP = DYN_COMP_LOSVD_MAP(bounds)

image = LOSVD_MAP[0]
image[0][0] = -1000

pc = axs[2, 2].imshow(image,cmap=cmap,transform = tr + axs[2,2].transData,vmin = Vel_min, vmax = Vel_max)
axs[2, 2].yaxis.set_major_locator(tcr.NullLocator())
axs[2, 2].xaxis.set_major_locator(tcr.NullLocator())

axs[2, 2].yaxis.set_major_formatter(tcr.NullFormatter())
axs[2, 2].xaxis.set_major_formatter(tcr.NullFormatter())
axs[2, 2].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[2, 2].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[2, 2].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 2].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 2].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 2].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 2].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[2, 2].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

image = LOSVD_MAP[1]
image[0][0] = -1000

pc = axs[2, 3].imshow(image,cmap=cmap,transform = tr + axs[2,3].transData,vmin = Sigma_min, vmax = Sigma_max)

axs[2, 3].yaxis.set_major_formatter(tcr.NullFormatter())
axs[2, 3].xaxis.set_major_formatter(tcr.NullFormatter())
axs[2, 3].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[2, 3].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[2, 3].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 3].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[2, 3].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 3].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2, 3].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[2, 3].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

######### Spherical component ###########

bounds = [-0.2,0.2]

image = orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(9,14) for radius in range(0,21)],0)
#image[(image < 1)] = None

pc = axs[3, 0].imshow(image,cmap=cmap_2,norm=norm)
axs[3, 0].yaxis.set_major_locator(tcr.NullLocator())
axs[3, 0].xaxis.set_major_locator(tcr.NullLocator())

axs[3, 0].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
axs[3, 0].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[3, 0].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[3, 0].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[3, 0].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 0].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 0].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 0].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 0].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[3, 0].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

axs[3, 0].annotate(r"$M_{c}/M_{total} =$" +  f"{int(np.sum(image) * 100)}" + "%",(1,30), fontsize=4,path_effects=[pe.withStroke(linewidth=1, foreground="white")])

axs[3, 0].annotate(f"{bounds[1]+0.05} > $\lambda_z$ > {bounds[0]-0.05}",(18,30), fontsize=3,path_effects=[pe.withStroke(linewidth=1, foreground="white")])

image = orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(9,14) for radius in range(0,21)],incl)
#image[(image < 1)] = None

pc = axs[3, 1].imshow(image,cmap=cmap_2,norm=norm)
axs[3, 1].yaxis.set_major_locator(tcr.NullLocator())
axs[3, 1].xaxis.set_major_locator(tcr.NullLocator())

axs[3, 1].yaxis.set_major_formatter(tcr.NullFormatter())
axs[3, 1].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[3, 1].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[3, 1].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[3, 1].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 1].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 1].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 1].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 1].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[3, 1].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

LOSVD_MAP = DYN_COMP_LOSVD_MAP(bounds)

image = LOSVD_MAP[0]
image[0][0] = -1000

pc = axs[3, 2].imshow(image,cmap=cmap,transform = tr + axs[3,2].transData,vmin = Vel_min, vmax = Vel_max)
axs[3, 2].yaxis.set_major_locator(tcr.NullLocator())
axs[3, 2].xaxis.set_major_locator(tcr.NullLocator())

axs[3, 2].yaxis.set_major_formatter(tcr.NullFormatter())
axs[3, 2].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[3, 2].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[3, 2].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[3, 2].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 2].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 2].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 2].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 2].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[3, 2].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

image = LOSVD_MAP[1]
image[0][0] = -1000

pc = axs[3, 3].imshow(image,cmap=cmap,transform = tr + axs[3,3].transData,vmin = Sigma_min, vmax = Sigma_max)
axs[3, 3].yaxis.set_major_locator(tcr.NullLocator())
axs[3, 3].xaxis.set_major_locator(tcr.NullLocator())

axs[3, 3].yaxis.set_major_formatter(tcr.NullFormatter())
axs[3, 3].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[3, 3].yaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))
axs[3, 3].xaxis.set_major_locator(tcr.IndexLocator(12,offset = -7.5))

axs[3, 3].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 3].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[3, 3].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 3].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[3, 3].tick_params(direction = "in",which = "major",width = 0.3,length = 2,labelsize =4)
axs[3, 3].tick_params(direction = "in",which = "minor",width = 0.15,length = 1)

fig.savefig("plots.pdf",format = "pdf")
