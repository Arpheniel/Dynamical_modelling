import numpy as np, agama as _agama
import matplotlib.pyplot as plt
import scipy
import sys
from time import sleep
from tqdm import tqdm
from decimal import Decimal
np.set_printoptions(threshold=sys.maxsize)
plt.get_cmap('inferno')
scale = 42
N = 16000  # number of orbits
M_total = 4.5e9  # in units of M_sun
ML = 15 #M/L
distance = 117490000# in [pc]
Filename = "M1e+07_O20_Rh150_Vh190_i42_a0_N16000_R1.00_GH_DensitySphHarm.npz"

archive = np.load(Filename, allow_pickle=True, encoding='latin1')

weights = archive["weights"]
lambda_z_list = archive["MOD_Lambda_z"]
Rmean_list = archive["Rmean"]
matrices = [archive["matrices1"],archive["matrices2"]]
cosincl = archive["cosincl"]

sorted_orbs = []

for i in range(20):
    sorted_orbs.append([])
    for j in range(40):
        sorted_orbs[i].append([])


for orbi in range(N):
    Rmean = Rmean_list[orbi]
    lambda_z = lambda_z_list[orbi]

    a = 1/3
    b = 10
    # i / a < r < (i + 1) / a
    # i / b + min(lambda_z_list) < lambda_z < (i + 1) / b + min(lambda_z_list)
    if int(Rmean * a) < 20 and int((lambda_z - min(lambda_z_list)) * b) < 40:
        sorted_orbs[int(Rmean * a)][int((lambda_z - min(lambda_z_list)) * b)].append(orbi)



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


def orbit_map(orb_groups):
    #width, height = 800, 800
    image = np.zeros((scale,scale)) 
    pixel_size = 1
    for orb_group in orb_groups:
        for orbi in orb_group:
            weight = weights[0][orbi]
            for str_point in open(f"orbits/orbit_{orbi}.txt"):
                #if weight >= 0.01: continue

                point = [float(i) for i in str_point.split(" ")]
                
                position = [point[0], point[1], point[2]]
                position = position @ rotation_matrix_x(np.radians(20))
                position = position @ rotation_matrix_y(np.radians(60))
                position = position @ rotation_matrix_z(np.radians(0))
                
                x_ = int(position[0] * scale*pixel_size*0.01) + int(scale/2*pixel_size)
                y_ = int(position[2] * scale*pixel_size*0.01) + int(scale/2*pixel_size)
                try:
                    image[y_//pixel_size][x_//pixel_size] += weight/1000#/(distance * np.pi / 648000)*(ML*M_total)
                except:
                    pass
    #fig,ax = plt.subplots()
    #ax.imshow(image)
    #fig.savefig("orb_map")
    return image


#orbit_map([sorted_orbs[radius][lambda_z] for lambda_z in range(0,40) for radius in range(0,10)])

from astropy.io import fits

file = fits.open("results_8992-3704_vorb020_md19_ad-1_nmom4.fits")

file.info()
binnum = file[5].data
age = file[22].data
age_err = file[23].data
met = file[24].data
met_err = file[25].data

bin_age = np.zeros(np.max(binnum)+1)
for x in range(0,42):
    for y in range(0,42):
        if binnum[x][y] >= 0:
            bin_age[binnum[x][y]] = met[x][y]

weight_matrix = []
dyn_comp_num = 0

weights_map = []
for i, lambda_z in zip(tqdm(range(0,40)), range(0,40)):
    for radius in range(0,20):
        if sorted_orbs[radius][lambda_z] == []: continue
        #print(sorted_orbs[radius][lambda_z])
        weights_map.append(orbit_map([sorted_orbs[radius][lambda_z]]))
        bin_orb_weight = np.zeros(np.max(binnum)+1)
        for x in range(0,42):
            for y in range(0,42):
                bin_orb_weight[binnum[x][y]] += weights_map[dyn_comp_num][x][y]
        dyn_comp_num += 1
        weight_matrix.append(bin_orb_weight)
    sleep(0.01)
print("##########")
print(np.sum(weight_matrix))
print(dyn_comp_num)
col_ind = []
row_ind = []

weight_array = []

for bin_num in range(0,np.max(binnum)+1):
    for din_comp in range(0,dyn_comp_num):
        col_ind.append(din_comp)
        row_ind.append(bin_num)
        norm_weight = weight_matrix[din_comp][bin_num]/np.sum(np.transpose(weight_matrix)[bin_num])
        if norm_weight > 0:
            weight_array.append(norm_weight)
        else:
            weight_array.append(0)
print(len(weight_array),len(col_ind),len(row_ind))
print(np.max(binnum)+1)



from scipy.sparse import csr_matrix


weight_matrix = csr_matrix((weight_array, (row_ind,col_ind)), shape=(np.max(binnum)+1,dyn_comp_num)).toarray()

"""
image = np.zeros((54,54))
for x in range(0,54):
    for y in range(0,54):
        for n in range(0,736):
            if binnum[x][y] == n:
                image[x][y] = np.mean(weight_matrix[n][100])

fig,ax = plt.subplots()
ax.imshow(image)
fig.savefig("bin_age")
"""

x = scipy.optimize.lsq_linear(weight_matrix, bin_age, bounds=(-2, 2), method='bvls', tol=1e-30, lsq_solver=None, lsmr_tol=None, max_iter=None, verbose=2, lsmr_maxiter=None)

print("###############")
print(x.x)


model_age = x.x

i = 0
young = [[],[],[]]
medium1 = [[],[],[]]
medium2 = [[],[],[]]
old = [[],[],[]]

for Rmean in range(0,20):
    for lambda_z in range(0,40):
        if sorted_orbs[Rmean][lambda_z] != []:
            num_of_dyn_comp = i
            if model_age[i] < -1:
                young[0].append(num_of_dyn_comp)
                young[1].append(lambda_z )
                young[2].append(Rmean)
            if model_age[i] > -1 and model_age[i] < -0.5:
                medium1[0].append(num_of_dyn_comp)
                medium1[1].append(lambda_z )
                medium1[2].append(Rmean)
            if model_age[i] > -0.5 and model_age[i] < 0:
                medium2[0].append(num_of_dyn_comp)
                medium2[1].append(lambda_z )
                medium2[2].append(Rmean)
            if model_age[i] > 0.5:
                old[0].append(num_of_dyn_comp)
                old[1].append(lambda_z )
                old[2].append(Rmean)
            i += 1

cm = 1/2.54  # centimeters in inches

fig,axs = plt.subplots(1,4,figsize=(20*cm, 10*cm))
num_comp = 0
dyn_comp_orbits = []
str_age_bound = ["0","4","6","8","10","12","14"]
str_met_bound = ["-2","-1","-0.5","0","0.5","1","2"]


for age in [young,medium1,medium2,old]:
    #stat_cosincl = 0
    image = np.zeros((scale,scale))
    for dyn_comp in range(0,len(age[0])):
        #dyn_comp_orbits.append(sorted_orbs[age[2][dyn_comp]][age[1][dyn_comp]])
        image += weights_map[age[0][dyn_comp]]
        #num_orb = 0
        #for orbi in sorted_orbs[age[2][dyn_comp]][age[1][dyn_comp]]:
         #   stat_cosincl += cosincl[orbi] * weights[0][orbi]
    #print(stat_cosincl/
    pc = axs[num_comp].imshow(image*M_total*ML,norm="log", cmap='inferno')
    axs[num_comp].set_xlabel("0.5 arcsec")
    axs[num_comp].set_ylabel("0.5 arcsec")
    M_t = '%.2E' % Decimal(str(np.sum(image)*M_total*ML))
    axs[num_comp].set_title(f'M_total = {M_t}', fontsize=8)
    axs[num_comp].text(6, -6, f'{str_met_bound[num_comp]} < Z < {str_met_bound[num_comp+1]}', fontsize=8)
    num_comp += 1
    dyn_comp_orbits = []


fig.colorbar(pc, ax=axs, orientation='horizontal', cmap='inferno',label ="M, M_sun")
fig.savefig("comp_weight_met" + str(num_comp))




"""
mean_age = np.zeros((42,42))
for x in range(0,42):
    for y in range(0,42):
        bin_weight_age = np.zeros(dyn_comp_num)
        for dyn_comp_i in range(0,dyn_comp_num):
            bin_weight_age[dyn_comp_i] = model_age[dyn_comp_i]*weight_matrix[binnum[x][y]][dyn_comp_i]
        mean_age[x][y] = np.mean(bin_weight_age)

fig,ax = plt.subplots()
ax.imshow(mean_age)
fig.savefig("model_met_map")

lambda_z_plot = []
Rmean_plot = []

for lambda_z in range(0,40):
    for Rmean in range(0,20):
        if sorted_orbs[Rmean][lambda_z] != []:
            lambda_z_plot.append(lambda_z/20 - 1)
            Rmean_plot.append(Rmean)
            print(sorted_orbs[Rmean][lambda_z])

fig,ax = plt.subplots()
pc = ax.scatter(lambda_z_plot,Rmean_plot, c = model_age, s=50, cmap='plasma')

fig.colorbar(pc, ax=ax, extend='both',label='age')
#fig.colorbar(label='age')
#ax.colorbar(x.x)
fig.savefig("model_met_plot")
"""
