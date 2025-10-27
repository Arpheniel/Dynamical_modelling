import numpy as np, agama as _agama
import matplotlib.pyplot as plt
from matplotlib import ticker as tcr
from matplotlib import transforms
from astropy.io import fits

M_total = 1.5e9 # in units of [M_sun]
Max_rad = 15 # in [arcsec]
Filename = "M1e+07_O0_Rh100_Vh190_i27_a0_N5000_R1.00_GH_DensityCylindricalLinear.npz"
R_scale = 0.5296 # scale of kpc in AGAMA
Im_scale = 32 # scale of Manga image (1 pixel = 0.5 arcsec)
distance = 108643 # in [kpc]
gamma = 0 # angle of rotation of losvd map
model_index = 0

archive = np.load(Filename, allow_pickle=True, encoding='latin1')

LOSVD = archive['LOSVD']
Lambda_z = archive['MOD_Lambda_z']
circ2 = archive["circ2"] # norm angular momentum
Vel = archive["Velocity"]
Upsilon = archive["Upsilon"]
weights = archive["weights"][model_index]
Rmean = archive["Rmean"]* R_scale/(distance * np.pi / 648000) # mean radius of orbit in [arcsec]
cosincl = archive["cosincl"]
ic = archive["ic"]
inttime = archive["inttime"]
DYN_COMP_LOSVD = archive["DYN_COMP_LOSVD"][model_index]

Units = _agama.getUnits()

spectr = fits.open("manga-8254-1902-MAPS-VOR10-MILESHC-MASTARSSP.fits")

phot = fits.open("mosaic-00157902-LEDA2220522-z-CCD3-image.fits")['CCD3'].data

def surf_br(surf, M_sun):
    return 10 ** (-(surf - M_sun - 21.57) / 2.5)


def magnitudo(flx):
    return 22.5 - 2.5 * np.log10(flx)
    
phot = magnitudo(phot)

phot = surf_br(phot, 4)

print(phot)

bin_scheme = np.loadtxt("bins_LEDA2220522.txt")



kinem_map = spectr['STELLAR_SIGMA'].data
kinem_map [(kinem_map  == 0)] = None

tr = transforms.Affine2D().rotate_deg_around(Im_scale/2,Im_scale/2,gamma)

def tickers_X_formatter(x, pos):
    #print(x,pos)
    return f'{abs((x - Im_scale/2)/2):.0f}'
def tickers_Y_formatter(y, pos):
    #print(y,pos)
    return f'{abs((y - Im_scale/2)/2):.0f}' #Y is inverted

cm = 1/2.54
fig, ax = plt.subplots(figsize=(14*cm, 11.5*cm))

cmap = plt.cm.nipy_spectral
cmap.set_bad(color='White') # Set 'bad' values (NaNs) to red

pc = ax.imshow(kinem_map,cmap = cmap,transform = tr + ax.transData,origin='lower',vmin = 0, vmax = 100)

ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
ax.yaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))
ax.xaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))

ax.yaxis.set_minor_formatter(tcr.NullFormatter())
ax.xaxis.set_minor_formatter(tcr.NullFormatter())
ax.yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
ax.xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

#ax.annotate(r"$M_{c}/M_{total} = 12$" + "%",(0.2,0.5), fontsize=15)
#ax.annotate(r"$M_{total} = 6.75 \cdot 10^{10} M_{\odot}$",(0,0.5),family='cursive', fontsize=15)

ax.annotate(r"$\sigma_{0}[km/s]$",(0,29),fontsize=15,style = "italic")
#ax.annotate(r"$V_{0}[km/s]$",(0,29),fontsize=15,style = "italic")
#ax.annotate(f"{Lambda_z_bound[1]} > $\lambda_z$ > {Lambda_z_bound[0]}",(25,39), fontsize=15)

ax.set_xlabel("arcsec")
ax.set_ylabel("arcsec")
#plt.title("PGC 35706" )
ax.tick_params(direction = "in",which = "both",labelsize = 8)
fig.colorbar(pc,ax=ax)
fig.savefig("obs_kinem_map.pdf",format = "pdf")



Lambda_z_bounds = [-1,1]

def DYN_COMP_LOSVD_MAP(bounds,GH_moment):

    LOSVD_plt = np.zeros((292,47))
    
    for r in range(21):
        for l_z in range(int((1 + Lambda_z_bounds[0])*10.5),int((1 + Lambda_z_bounds[1])*10.5)): # 0 = -1, 11 = 0 , 21 = 1
            LOSVD_plt +=  DYN_COMP_LOSVD[l_z][r]

    GH_moments = _agama.ghMoments(matrix=LOSVD_plt * 7**-0.5,gridv=np.linspace(-250, 250, 46) * 7**0.5, degree=2, ghorder=6)[:,(1,2,6,7,8,9)]

    kinem_map = np.full((Im_scale,Im_scale),None,dtype=float)

    for bin_i in range(len(bin_scheme)):
        kinem_map[int(bin_scheme[bin_i][0]+Im_scale/2)][int(bin_scheme[bin_i][1]+Im_scale/2)] = GH_moments[int(bin_scheme[bin_i][2])][GH_moment]
    
    return kinem_map

tr = transforms.Affine2D().rotate_deg_around(Im_scale/2,Im_scale/2,gamma)

def tickers_X_formatter(x, pos):
    return f'{abs((x - Im_scale/2)/2):.0f}'
def tickers_Y_formatter(y, pos): #for kinem
    return f'{abs((y - Im_scale/2)/2):.0f}'
    
def tickers_Y_formatter_2(y, pos): #for phot
    return f'{abs((y - Im_scale/2)+4):.0f}'

#cm = 1/2.54
#fig, ax = plt.subplots(figsize=(14*cm, 11.5*cm))

cmap = plt.cm.nipy_spectral
cmap.set_bad(color='White') # Set 'bad' values (NaNs) to white
cmap_2 = plt.cm.inferno_r

#kinem_map = kinem_map - spectr['STELLAR_VEL'].data

Vel_map_obs = spectr['STELLAR_VEL'].data
Vel_map_obs[(Vel_map_obs == 0)] = None

fig, axs = plt.subplots(3, 3, sharex = False, sharey = False)

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=-0.58)

### Vel obs vs Vel model ###

Vel_min = -45
Vel_max = 45

Sig_min = 0
Sig_max = 100

pc = axs[1,0].imshow(Vel_map_obs,cmap = cmap,transform = tr + axs[1,0].transData,origin='lower',vmin = Vel_min, vmax = Vel_max)

axs[1, 0].yaxis.set_major_locator(tcr.NullLocator())
axs[1, 0].xaxis.set_major_locator(tcr.NullLocator())

#axs[1, 0].yaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 0].xaxis.set_major_formatter(tcr.NullFormatter())

axs[1,0].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
#axs[1,0].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[1,0].yaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))
axs[1,0].xaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))

axs[1,0].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[1,0].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[1,0].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1,0].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

axs[1,0].set_ylabel("arcsec")

#ax.annotate(r"$M_{c}/M_{total} = 12$" + "%",(0.2,0.5), fontsize=15)
#ax.annotate(r"$M_{total} = 6.75 \cdot 10^{10} M_{\odot}$",(0,0.5),family='cursive', fontsize=15)

axs[1,0].annotate(r"$V_{0}[km/s]$",(2,28),fontsize=7,style = "italic")
#axs[1,0].annotate(f"{Lambda_z_bounds[1]} > $\lambda_z$ > {Lambda_z_bounds[0]}",(25,39), fontsize=15)

#axs[1,0].set_xlabel("arcsec")
#axs[1,0].set_ylabel("arcsec")
#plt.title("PGC_35706" )
axs[1,0].tick_params(direction = "in",which = "both",labelsize = 6)
#fig.colorbar(pc,ax=ax)

kinem_map = DYN_COMP_LOSVD_MAP(Lambda_z_bounds ,0)

pc = axs[1,1].imshow(kinem_map,cmap = cmap,transform = tr + axs[1,1].transData,origin='lower',vmin = Vel_min, vmax = Vel_max)

axs[1, 1].yaxis.set_major_locator(tcr.NullLocator())
axs[1, 1].xaxis.set_major_locator(tcr.NullLocator())

axs[1, 1].yaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 1].xaxis.set_major_formatter(tcr.NullFormatter())

#axs[1,1].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
#axs[1,1].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[1,1].yaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))
axs[1,1].xaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))

axs[1,1].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[1,1].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[1,1].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1,1].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

axs[1,1].tick_params(direction = "in",which = "both",labelsize = 6)



kinem_map = kinem_map - spectr['STELLAR_VEL'].data

pc = axs[1,2].imshow(kinem_map,cmap = cmap,transform = tr + axs[1,2].transData,origin='lower',vmin = Vel_min, vmax = Vel_max)

divider = make_axes_locatable(axs[1,2])
ax_cb = divider.append_axes("right", size="5%", pad=-0.06)

fig.colorbar(pc, ax=axs[1, 2],cax = ax_cb,orientation = 'vertical',location = 'right',ticks = [Vel_min,Vel_max],format = tcr.NullFormatter())

axs[1,2].yaxis.set_label_position("right")
axs[1,2].set_ylabel(f'{Vel_min}                 {Vel_max} ',loc = 'top',labelpad = 1)



axs[1, 2].yaxis.set_major_locator(tcr.NullLocator())
axs[1, 2].xaxis.set_major_locator(tcr.NullLocator())

axs[1, 2].yaxis.set_major_formatter(tcr.NullFormatter())
axs[1, 2].xaxis.set_major_formatter(tcr.NullFormatter())

#axs[1,2].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
#axs[1,2].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[1,2].yaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))
axs[1,2].xaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))

axs[1,2].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[1,2].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[1,2].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[1,2].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

axs[1,2].tick_params(direction = "in",which = "both",labelsize = 6)



### Phot obs vs Phot model ###


phot_obs = np.loadtxt('converted_fits_LEDA_2220522_z.txt')
phot_model = np.loadtxt('model_LEDA_2220522_z.txt')
phot_res = np.abs(phot_obs - phot_model)
#phot_res = phot_obs - phot_model

Phot_max = np.max([np.max(phot_obs),np.max(phot_model),np.max(phot_res)])
Phot_min = np.min([np.min(phot_obs),np.min(phot_model),np.min(phot_res)])

pc = axs[0,0].imshow(phot_obs,cmap = cmap_2,transform = tr + axs[0,0].transData,origin='lower',vmin = Phot_min, vmax = Phot_max)
print(phot_obs.shape)
axs[0, 0].yaxis.set_major_locator(tcr.NullLocator())
axs[0, 0].xaxis.set_major_locator(tcr.NullLocator())

#axs[0, 0].yaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 0].xaxis.set_major_formatter(tcr.NullFormatter())

axs[0,0].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter_2))
#axs[0,0].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[0,0].yaxis.set_major_locator(tcr.IndexLocator(5,offset = 2.5))
axs[0,0].xaxis.set_major_locator(tcr.IndexLocator(5,offset = 2.5))

axs[0,0].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[0,0].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[0,0].yaxis.set_minor_locator(tcr.IndexLocator(1,offset = 0.5))
axs[0,0].xaxis.set_minor_locator(tcr.IndexLocator(1,offset = 0.5))

axs[0,0].tick_params(direction = "in",which = "both",labelsize = 6)

axs[0,0].annotate(r"$Log_{10}(Flux)[counts/pixel]$",(1,23),fontsize=5,style = "italic")

axs[0,0].set_title("obs ")



pc = axs[0,1].imshow(phot_model,cmap = cmap_2,transform = tr + axs[0,1].transData,origin='lower',vmin = Phot_min, vmax = Phot_max)

axs[0, 1].yaxis.set_major_locator(tcr.NullLocator())
axs[0, 1].xaxis.set_major_locator(tcr.NullLocator())

axs[0, 1].yaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 1].xaxis.set_major_formatter(tcr.NullFormatter())

#axs[0,1].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
#axs[0,1].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[0,1].yaxis.set_major_locator(tcr.IndexLocator(5,offset = 2.5))
axs[0,1].xaxis.set_major_locator(tcr.IndexLocator(5,offset = 2.5))

axs[0,1].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[0,1].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[0,1].yaxis.set_minor_locator(tcr.IndexLocator(1,offset = 0.5))
axs[0,1].xaxis.set_minor_locator(tcr.IndexLocator(1,offset = 0.5))

axs[0,1].tick_params(direction = "in",which = "both",labelsize = 6)

axs[0,1].set_title("model ")



pc = axs[0,2].imshow(phot_res,cmap = cmap_2,transform = tr + axs[0,2].transData,origin='lower',vmin = Phot_min, vmax = Phot_max)

divider = make_axes_locatable(axs[0,2])
ax_cb = divider.append_axes("right", size="5%", pad=-0.06)

fig.colorbar(pc, ax=axs[0, 2],cax = ax_cb,orientation = 'vertical',location = 'right',ticks = [Phot_min,Phot_max],format = tcr.NullFormatter())

axs[0,2].yaxis.set_label_position("right")
axs[0,2].set_ylabel(f'{round(Phot_min)}                  {round(Phot_max)} ',loc = 'top',labelpad = 1)

axs[0, 2].yaxis.set_major_locator(tcr.NullLocator())
axs[0, 2].xaxis.set_major_locator(tcr.NullLocator())

axs[0, 2].yaxis.set_major_formatter(tcr.NullFormatter())
axs[0, 2].xaxis.set_major_formatter(tcr.NullFormatter())

#axs[0,2].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
#axs[0,2].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[0,2].yaxis.set_major_locator(tcr.IndexLocator(5,offset = 2.5))
axs[0,2].xaxis.set_major_locator(tcr.IndexLocator(5,offset = 2.5))

axs[0,2].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[0,2].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[0,2].yaxis.set_minor_locator(tcr.IndexLocator(1,offset = 0.5))
axs[0,2].xaxis.set_minor_locator(tcr.IndexLocator(1,offset = 0.5))

axs[0,2].tick_params(direction = "in",which = "both",labelsize = 6)

axs[0, 2].set_title("|obs - model|")



### Sigma obs vs Sigma model ###



Sig_map_obs = spectr['STELLAR_SIGMA'].data
Sig_map_obs[(Sig_map_obs == 0)] = None

pc = axs[2,0].imshow(Sig_map_obs,cmap = cmap,transform = tr + axs[2,0].transData,origin='lower',vmin = Sig_min, vmax = Sig_max)

axs[2,0].yaxis.set_major_locator(tcr.NullLocator())
axs[2,0].xaxis.set_major_locator(tcr.NullLocator())

#axs[2,0].yaxis.set_major_formatter(tcr.NullFormatter())
#axs[2,0].xaxis.set_major_formatter(tcr.NullFormatter())

axs[2,0].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
axs[2,0].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[2,0].yaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))
axs[2,0].xaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))

axs[2,0].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[2,0].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[2,0].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2,0].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

axs[2,0].tick_params(direction = "in",which = "both",labelsize = 6)

axs[2,0].annotate(r"$\sigma_{0}[km/s]$",(2,28),fontsize=7,style = "italic")



kinem_map = DYN_COMP_LOSVD_MAP(Lambda_z_bounds ,1)


pc = axs[2,1].imshow(kinem_map,cmap = cmap,transform = tr + axs[2,1].transData,origin='lower',vmin = Sig_min, vmax = Sig_max)

axs[2,1].yaxis.set_major_locator(tcr.NullLocator())
axs[2,1].xaxis.set_major_locator(tcr.NullLocator())

axs[2,1].yaxis.set_major_formatter(tcr.NullFormatter())
#axs[2,1].xaxis.set_major_formatter(tcr.NullFormatter())

#axs[2,1].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
axs[2,1].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[2,1].yaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))
axs[2,1].xaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))

axs[2,1].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[2,1].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[2,1].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2,1].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

axs[2,1].tick_params(direction = "in",which = "both",labelsize = 6)

axs[2,1].set_xlabel("arcsec")


kinem_map = abs(kinem_map - spectr['STELLAR_SIGMA'].data)


pc = axs[2,2].imshow(kinem_map,cmap = cmap,transform = tr + axs[2,2].transData,origin='lower',vmin = Sig_min, vmax = Sig_max)

divider = make_axes_locatable(axs[2,2])
ax_cb = divider.append_axes("right", size="5%", pad=-0.06)

fig.colorbar(pc, ax=axs[2, 2],cax = ax_cb,orientation = 'vertical',location = 'right',ticks = [Sig_min,Sig_max],format = tcr.NullFormatter())

axs[2,2].yaxis.set_label_position("right")
axs[2,2].set_ylabel(f'{round(Sig_min)}                  {round(Sig_max)} ',loc = 'top',labelpad = 1)

axs[2,2].yaxis.set_major_locator(tcr.NullLocator())
axs[2,2].xaxis.set_major_locator(tcr.NullLocator())

axs[2,2].yaxis.set_major_formatter(tcr.NullFormatter())
#axs[2,2].xaxis.set_major_formatter(tcr.NullFormatter())

#axs[2,2].yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
axs[2,2].xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
axs[2,2].yaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))
axs[2,2].xaxis.set_major_locator(tcr.IndexLocator(10,offset = 6.5))

axs[2,2].yaxis.set_minor_formatter(tcr.NullFormatter())
axs[2,2].xaxis.set_minor_formatter(tcr.NullFormatter())
axs[2,2].yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))
axs[2,2].xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 0.5))

axs[2,2].tick_params(direction = "in",which = "both",labelsize = 6)

fig.savefig("kinem_map.pdf",format = "pdf")




orbits=[]
mass_c_p=[]
mass_c_n=[]
mass_sigma=[]
mass_polar_ring=[]
M_p=0
M_n=0
M_s=0
M_p_r=0

N = len(circ2)

for i in range(0,N):
    if weights[i]>=0.00001:
        orbits.append([circ2[i],weights[i] * M_total*7/1e8,Rmean[i] * R_scale,cosincl[i], Lambda_z[i], Vel[i]])

steps = 21
rad = np.linspace(0,Max_rad,steps)
zone = Max_rad/steps/2

for r in rad:
    for i in range(0,len(orbits)):
        if orbits[i][4] > 0.5 and orbits[i][2] > (r-zone) and orbits[i][2] < (r+zone):
            M_p += orbits[i][1] /10
        if orbits[i][4] < -0.5 and orbits[i][2] > (r-zone) and orbits[i][2] < (r+zone):
            M_n += orbits[i][1] /10
        if orbits[i][2] > (r-zone) and orbits[i][2] < (r+zone):
            M_s += orbits[i][1] /10
        if orbits[i][2] > (r-zone) and orbits[i][2] < (r+zone) and orbits[i][4] > -0.5 and orbits[i][4] < 0.5:
            M_p_r += orbits[i][1] /10
        
    print("#################",M_s,M_p,M_n,M_p_r)
    mass_sigma.append(M_s)
    mass_c_p.append(M_p)
    mass_c_n.append(M_n)
    mass_polar_ring.append(M_p_r)
    M_p = 0
    M_n = 0
    M_s = 0
    M_p_r = 0

fig, ax = plt.subplots()

ax.plot(rad,mass_c_p,label="mass of co-rotating component",marker='o', markersize = 4, linestyle="--",c = "blue")

ax.plot(rad,mass_c_n,label="mass of counter-rotating component",marker='o', markersize = 4 ,linestyle="--",c = "red")

ax.plot(rad,mass_sigma,label="mass of all stars",marker='o', markersize = 4 ,linestyle="--",c = "black")

ax.plot(rad,mass_polar_ring,label="mass of polar rings",marker='o', markersize = 4 ,linestyle="--",c = "grey")

ax.set_xlabel(r"$R_{mean},arcsec$")
ax.set_ylabel(r"Mass,$10^9$ $M_{\odot}$")
ax.set_title(r"$M_{Total}=9 \cdot 10^{9}$ $M_{\odot}$" )
ax.tick_params(direction = "in",which = "both",labelsize = 6)

fig.savefig("mass_dest.pdf",format = "pdf")


def orb_weights_dist(X,Y,orbits_index_1,orbits_index_2,steps,X_label,Y_label,Plot_title,plot_Name): 
# X - linspace of X , Y - linspace of Y, orbits_index_1 - index of X in orbits matrix, orbits_index_2 - index of Y in orbits matrix, steps - linspace steps
    
    Max_X , Max_Y = np.max(X) , np.max(Y)
    image_orb_weights_dist = np.zeros((steps,steps))
    
    for x in X:
        for y in Y:
            for orbi in range(0,len(orbits)):
                if orbits[orbi][orbits_index_1] > x - (Max_X/steps/2) and orbits[orbi][orbits_index_1] < x + (Max_X/steps/2) and orbits[orbi][orbits_index_2] > y - (Max_Y/steps/2) and orbits[orbi][orbits_index_2] < y + (Max_Y/steps/2):
                    image_orb_weights_dist[steps - 1 - int(round(((y+np.abs(np.min(Y)))/(Max_Y+np.abs(np.min(Y))))*(steps-1)))][int(round(((x+np.abs(np.min(X)))/(Max_X+np.abs(np.min(X))))*(steps-1)))] += orbits[orbi][1]

    def tickers_X_formatter(x, pos):
        return f'{x/((steps-1)/(np.abs(np.min(X))+Max_X))+np.min(X):.0f}'

    def tickers_Y_formatter(y, pos):
        return f'{Max_Y - np.abs(np.min(Y)) -(y/((steps-1)/(np.abs(np.min(Y))+Max_Y))+np.min(Y)):.1f}' #Y is inverted

    fig, ax = plt.subplots()
    pc = ax.imshow(image_orb_weights_dist,cmap = "inferno")
    ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
    ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
    ax.xaxis.set_minor_formatter(tcr.NullFormatter())
    ax.yaxis.set_major_locator(tcr.IndexLocator(2,offset = 0.4))
    ax.xaxis.set_minor_locator(tcr.IndexLocator(1,offset = 0.5))
    ax.set_xlabel(X_label)
    ax.set_ylabel(Y_label)
    fig.colorbar(pc,ax=ax,label=r'$Mass, 10^8$ $M_{\odot}$')
    fig.savefig(plot_Name + ".pdf",format = "pdf")

steps = 21

rad = np.linspace(Max_rad,0,steps)
cos_incl = np.linspace(1,-1,steps)

orb_weights_dist(rad,cos_incl,2,3,steps,"R_mean,kpc","cos_incl","","Rad_vs_incl")

norm_L = np.linspace(1,0,steps)

orb_weights_dist(rad,norm_L,2,0,steps,"R_mean,kpc","norm_L","","Rad_vs_norm_L")

plt_Lambda_z = np.linspace(-1,1,steps)

orb_weights_dist(rad,plt_Lambda_z,2,4,steps,r"$R_{mean},arcsec$",r"$circularity,\lambda_{z}$","","Rad_vs_Lambda_z")
