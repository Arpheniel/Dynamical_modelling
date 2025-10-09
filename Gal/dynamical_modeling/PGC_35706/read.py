import numpy as np, agama as _agama
import matplotlib.pyplot as plt
from matplotlib import ticker as tcr
from matplotlib import transforms

M_total = 4.5e9 # in units of [M_sun]
Max_rad = 20 # in [kpc]
Filename = "M1e+07_O0_Rh169_Vh185_i42_a0_N6000_R1.00_GH_DensitySphHarm.npz"
R_scale = 0.5696 # scale of kpc in AGAMA
Im_scale = 42
distance = 117490

archive = np.load(Filename, allow_pickle=True, encoding='latin1')

LOSVD = archive['LOSVD']
Lambda_z = archive['MOD_Lambda_z']
circ2 = archive["circ2"]
Vel = archive["Velocity"]
Upsilon = archive["Upsilon"]
weights = archive["weights"]
Rmean = archive["Rmean"]/(distance * np.pi / 648000)
cosincl = archive["cosincl"]
ic = archive["ic"]
inttime = archive["inttime"]
DYN_COMP_LOSVD = archive["DYN_COMP_LOSVD"]

Units = _agama.getUnits()

bin_scheme = np.loadtxt("bins_PGC35706.txt")

DYN_COMP_LOSVD = DYN_COMP_LOSVD[0]

Lambda_z_bound = [0.5,1]

def DYN_COMP_LOSVD_MAP(bounds):

    LOSVD_plt = np.zeros((695,47))
    
    for r in range(21):
        for l_z in range(int((1 + bounds[0])*10.5),int((1 + bounds[1])*10.5)): # 0 = -1, 11 = 0 , 21 = 1
            LOSVD_plt +=  DYN_COMP_LOSVD[l_z][r]
    return LOSVD_plt
    
LOSVD_plt = DYN_COMP_LOSVD_MAP(Lambda_z_bound)

#LOSVD = LOSVD[0][0]

GH_moments = _agama.ghMoments(matrix=LOSVD_plt * 15**-0.5,gridv=np.linspace(-250, 250, 46) * 15**0.5, degree=2, ghorder=6)[:,(1,2,6,7,8,9)]

kinem_map = np.full((Im_scale,Im_scale),None,dtype=float)



for bin_i in range(len(bin_scheme)):
    kinem_map[int(bin_scheme[bin_i][0]+Im_scale/2)][int(bin_scheme[bin_i][1]+Im_scale/2)] = GH_moments[int(bin_scheme[bin_i][2])][0]

gamma = 0

tr = transforms.Affine2D().rotate_deg_around(Im_scale/2,Im_scale/2,gamma)

def tickers_X_formatter(x, pos):
    #print(x,pos)
    return f'{(x - Im_scale/2)/2:.1f}'
def tickers_Y_formatter(y, pos):
    #print(y,pos)
    return f'{(y - Im_scale/2)/2:.1f}' #Y is inverted

cm = 1/2.54
fig, ax = plt.subplots(figsize=(14*cm, 11.5*cm))

cmap = plt.cm.nipy_spectral
cmap.set_bad(color='White') # Set 'bad' values (NaNs) to red

pc = ax.imshow(kinem_map,cmap = cmap,transform = tr + ax.transData,origin='lower',vmin = -200, vmax = 200)

ax.yaxis.set_major_formatter(tcr.FuncFormatter(tickers_Y_formatter))
ax.xaxis.set_major_formatter(tcr.FuncFormatter(tickers_X_formatter))
ax.yaxis.set_major_locator(tcr.IndexLocator(10,offset = 1.5))
ax.xaxis.set_major_locator(tcr.IndexLocator(10,offset = 1.5))

ax.yaxis.set_minor_formatter(tcr.NullFormatter())
ax.xaxis.set_minor_formatter(tcr.NullFormatter())
ax.yaxis.set_minor_locator(tcr.IndexLocator(2,offset = 1.5))
ax.xaxis.set_minor_locator(tcr.IndexLocator(2,offset = 1.5))

ax.annotate(r"$M_{c}/M_{total} = 12$" + "%",(0.2,0.5), fontsize=15)
#ax.annotate(r"$M_{total} = 6.75 \cdot 10^{10} M_{\odot}$",(0,0.5),family='cursive', fontsize=15)

ax.annotate(r"$V_{0}[km/s]$",(0,38),fontsize=15,style = "italic")
ax.annotate(f"{Lambda_z_bound[1]} > $\lambda_z$ > {Lambda_z_bound[0]}",(25,39), fontsize=15)

ax.set_xlabel("arcsec")
ax.set_ylabel("arcsec")
plt.title("PGC_35706" )
ax.tick_params(direction = "in",which = "both")
fig.colorbar(pc,ax=ax)
fig.savefig("kinem_map",format = "eps")


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
    if weights[0][i]>=0.00001:
        orbits.append([circ2[i],weights[0][i] * M_total*15/1e8,Rmean[i] * R_scale,cosincl[i], Lambda_z[i], Vel[i]])

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

matrix=[[],[],[],[],[],[],[]]

for i in range(0,len(orbits)):
        matrix[0].append(orbits[i][0])
        matrix[1].append(orbits[i][1])
        matrix[2].append(orbits[i][2])
        matrix[3].append(orbits[i][3])
        matrix[4].append(orbits[i][4])
        matrix[5].append(orbits[i][5])

fig, ax = plt.subplots()

ax.plot(rad,mass_c_p,label="mass of co-rotating component",marker='o', markersize = 4, linestyle="--",c = "blue")

ax.plot(rad,mass_c_n,label="mass of counter-rotating component",marker='o', markersize = 4 ,linestyle="--",c = "red")

ax.plot(rad,mass_sigma,label="mass of all stars",marker='o', markersize = 4 ,linestyle="--",c = "black")

ax.plot(rad,mass_polar_ring,label="mass of polar rings",marker='o', markersize = 4 ,linestyle="--",c = "grey")

#ax.legend()
ax.set_xlabel(r"$R_{mean},arcsec$")
ax.set_ylabel(r"Mass,$10^9$ $M_{\odot}$")
ax.set_title(r"$M_{Total}=6.75 \cdot 10^{10}$ $M_{\odot}$" )
#ax.set_ylim(0, 10)
#ax.set_yscale('log')
#ax.yaxis.set_major_locator(tcr.LogLocator(base = 10, numticks = 200))
ax.tick_params(direction = "in",which = "both")
#ax.xaxis.set_minor_formatter(tcr.NullFormatter())
#ax.xaxis.set_minor_locator(tcr.IndexLocator(0.5,offset = 0))
fig.savefig("mass_dest",format = "eps")


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
    fig.savefig(plot_Name,format = "eps")

steps = 21

rad = np.linspace(Max_rad,0,steps)
cos_incl = np.linspace(1,-1,steps)

orb_weights_dist(rad,cos_incl,2,3,steps,"R_mean,kpc","cos_incl","","Rad_vs_incl")

norm_L = np.linspace(1,0,steps)

orb_weights_dist(rad,norm_L,2,0,steps,"R_mean,kpc","norm_L","","Rad_vs_norm_L")

plt_Lambda_z = np.linspace(-1,1,steps)

orb_weights_dist(rad,plt_Lambda_z,2,4,steps,r"$R_{mean},arcsec$",r"$circularity,\lambda_{z}$","","Rad_vs_Lambda_z")

print(np.mean(Rmean),np.max(Rmean))
"""
matrix_disk = [[],[],[],[],[],[],[]]

for i in range(0,len(orbits)):
    if orbits[i][4] > 0.7:
        matrix_disk[0].append(orbits[i][0])
        matrix_disk[1].append(orbits[i][1])
        matrix_disk[2].append(orbits[i][2])
        matrix_disk[3].append(orbits[i][3])
        matrix_disk[4].append(orbits[i][4])
        matrix_disk[5].append(orbits[i][5])
        matrix_disk[6].append(orbits[i][6])
        
fig, ax = plt.subplots()
ax.scatter(matrix_disk[2],matrix_disk[6],20,c = "red")
ax.set_xlabel("R_mean,kpc")
ax.set_ylabel("Vel,km/s")
fig.savefig("Vel_vs_Rad")

fig, axs = plt.subplots(2, 2, layout='constrained')
pc = axs[0, 0].scatter(matrix[1], matrix[2], c=matrix[3], s=0.01 ,cmap='inferno',alpha = 1)
fig.colorbar(pc, ax=axs[0, 0], extend='both',label='cos(incl)')
#axs[0, 0].set_title('c=cos(incl)')
axs[0, 0].set_xlabel('weights,M_sun')
axs[0, 0].set_ylabel('R_mean,kpc')

pc = axs[1, 0].scatter(matrix[1], matrix[2], c=matrix[0], s=10 ,cmap='inferno',alpha = 0.2)
fig.colorbar(pc, ax=axs[1, 0], extend='both',label='circ2')
#axs[1, 0].set_title('c=circ2')
axs[1, 0].set_xlabel('weights,M_sun')
axs[1, 0].set_ylabel('R_mean,kpc')

pc = axs[0, 1].scatter(matrix[2], matrix[3], c=matrix[1], s=10 ,cmap='inferno',alpha = 0.2)
fig.colorbar(pc, ax=axs[0, 1], extend='both',label='weights,M_sun')
#axs[0, 1].set_title('c=weights,M_sun')
axs[0, 1].set_xlabel('R_mean,kpc')
axs[0, 1].set_ylabel('cos(incl)')

pc = axs[1, 1].scatter(matrix[1], matrix[0], c=matrix[3], s=10 ,cmap='inferno',alpha = 0.2)
fig.colorbar(pc, ax=axs[1, 1], extend='both',label='cos(incl)')
#axs[1, 1].set_title('c=cos(incl)')
axs[1, 1].set_xlabel('weights,M_sun')
axs[1, 1].set_ylabel('circ2')

fig.savefig("plots")

#data = _numpy.transpose(_numpy.array([circ2,weights[5],Rmean,cosincl]))
print(Units)
"""
"""
np.savetxt("/home/denis/Documents/agama/Gal/data.txt",sorted(orbits, key = lambda orbits: orbits[1]), fmt="%8.3f", header="circ2 weights Rmean Lz/l")
"""
