import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from os import path
import vorbin
import sys
import argparse
from vorbin.voronoi_2d_binning import voronoi_2d_binning

spectro = fits.open("NGC524_stell_cont_5.0_20.0_ppxf_results_emiles_4700_9000.fits")[1] # open MUSE
bins_file = fits.open("NGC524_stell_cont_5.0_20.0_binning.fits")

#LOF_vel = 7179 #line of sight velocity

bin_list = bins_file[1].data
bin_x = bins_file[2].data
bin_y = bins_file[3].data

x_scale = bins_file[0].shape[0]
y_scale = bins_file[0].shape[1]

x = []
y = []

for x_ in range(0,x_scale+1):
    for y_ in range(0,y_scale+1):
        x.append(x_)
        y.append(y_)
        
x = np.array(x)
y = np.array(y)

xBin = bin_x
yBin = bin_y

if x.size < 1e4:
    binNum = np.argmin((x[:, None] - xBin)**2 + (y[:, None] - yBin)**2, axis=1)
else:  # use for loop to reduce memory usage
    binNum = np.zeros(x.size, dtype=int)
    for j, (xj, yj) in enumerate(zip(x, y)):
        binNum[j] = np.argmin((xj - xBin)**2 + (yj - yBin)**2)
bin_id = np.transpose([x - x_scale/2,y- y_scale/2,binNum])

np.savetxt('bins_NGC_524.txt', bin_id,fmt='%8.4f %8.4f %7d',
           header='Binning scheme file\nX_pix  Y_pix  binid')



VEL_map = spectro.data["V"]
VEL_ERR_map = spectro.data["E_V"]
SIG_map = spectro.data["SIG"]
SIG_ERR_map = spectro.data["E_SIG"]
H3_map = spectro.data["H3"]
H3_ERR_map = spectro.data["E_H3"]
H4_map = spectro.data["H4"]
H4_ERR_map = spectro.data["E_H4"]
H5_map = spectro.data["H5"]
H5_ERR_map = spectro.data["E_H5"]
H6_map = spectro.data["H6"]
H6_ERR_map = spectro.data["E_H6"]



VEL_map_1 = []
SIG_map_1 = []
H3_map_1 = []
H4_map_1 = []
H5_map_1 = []
H6_maP_1 = []
VEL_ERR_map_1 = []
SIG_ERR_map_1 = []
H3_ERR_map_1 = []
H4_ERR_map_1 = []
H5_ERR_map_1 = []
H6_ERR_map_1 =[]

VEL_map_2 = []
SIG_map_2 = []
H3_map_2 = []
H4_map_2 = []
H5_map_2 = []
H6_maP_2 = []
VEL_ERR_map_2 = []
SIG_ERR_map_2 = []
H3_ERR_map_2 = []
H4_ERR_map_2 = []
H5_ERR_map_2 = []
H6_ERR_map_2 = []

for map_i in range(12402):
    VEL_map_1.append(VEL_map[map_i][0][0])
    SIG_map_1.append(SIG_map[map_i][0][0])
    H3_map_1.append(H3_map[map_i][0][0])
    H4_map_1.append(H4_map[map_i][0][0])
    H5_map_1.append(H5_map[map_i][0][0])
    H6_maP_1.append(H6_map[map_i][0][0])
    VEL_ERR_map_1.append(VEL_ERR_map[map_i][0][0])
    SIG_ERR_map_1.append(SIG_ERR_map[map_i][0][0])
    H3_ERR_map_1.append(H3_ERR_map[map_i][0][0])
    H4_ERR_map_1.append(H4_ERR_map[map_i][0][0])
    H5_ERR_map_1.append(H5_ERR_map[map_i][0][0])
    H6_ERR_map_1.append(H6_ERR_map[map_i][0][0])
    
    VEL_map_2.append(VEL_map[map_i][1][0])
    SIG_map_2.append(SIG_map[map_i][1][0])
    H3_map_2.append(H3_map[map_i][1][0])
    H4_map_2.append(H4_map[map_i][1][0])
    H5_map_2.append(H5_map[map_i][1][0])
    H6_maP_2.append(H6_map[map_i][1][0])
    VEL_ERR_map_2.append(VEL_ERR_map[map_i][1][0])
    SIG_ERR_map_2.append(SIG_ERR_map[map_i][1][0])
    H3_ERR_map_2.append(H3_ERR_map[map_i][1][0])
    H4_ERR_map_2.append(H4_ERR_map[map_i][1][0])
    H5_ERR_map_2.append(H5_ERR_map[map_i][1][0])
    H6_ERR_map_2.append(H6_ERR_map[map_i][1][0])

LOS_vel = np.mean(VEL_map) #line of sight velocity

#nump_data_1 = np.transpose([VEL_map_1 - LOS_vel,VEL_ERR_map_1,SIG_map_1,SIG_ERR_map_1,H3_map_1,H3_ERR_map_1,H4_map_1,H4_ERR_map_1])
nump_data_1 = np.transpose([VEL_map_1 - LOS_vel,VEL_ERR_map_1,SIG_map_1,SIG_ERR_map_1])

np.savetxt('kinem_gh_NGC_524.txt', nump_data_1, fmt='%8.3f',
           header='v        v_err    sigma    sigma_err')
