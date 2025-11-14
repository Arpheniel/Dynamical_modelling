import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from os import path
import vorbin
import sys
import argparse
from vorbin.voronoi_2d_binning import voronoi_2d_binning

spectro = fits.open("results_8992-3704_vorb020_md19_ad-1_nmom4.fits")  # open MANGA results of analysis

#LOF_vel = 7179 #line of sight velocity

spectro.info()
Bin_map = spectro[5].data
Vel_map = spectro[6].data
VEL_ERR_map = spectro[7].data
SIG_map = spectro[8].data
SIG_ERR_map = spectro[9].data
H3_map = spectro[10].data
H3_ERR_map = spectro[11].data
H4_map = spectro[12].data
H4_ERR_map = spectro[13].data

Len_bin = np.max(Bin_map) + 1

VEL_BIN = np.zeros(Len_bin)
VEL_ERR_BIN = np.zeros(Len_bin)
SIG_BIN = np.zeros(Len_bin)
SIG_ERR_BIN = np.zeros(Len_bin)
H3_BIN = np.zeros(Len_bin)
H3_ERR_BIN = np.zeros(Len_bin)
H4_BIN = np.zeros(Len_bin)
H4_ERR_BIN = np.zeros(Len_bin)

BINS = []

shape = Bin_map.shape[0]

for x in range(0,shape):
    for y in range(0,shape):
        for Bin_i in range(0,Len_bin):
            if Bin_map[x][y] == Bin_i:
                BINS.append([x-shape/2,y-shape/2,Bin_i])
                VEL_BIN[Bin_i] = Vel_map[x][y]
                SIG_BIN[Bin_i] = SIG_map[x][y]
                H3_BIN[Bin_i] = H3_map[x][y]
                H4_BIN[Bin_i] = H4_map[x][y]
                VEL_ERR_BIN[Bin_i] = VEL_ERR_map[x][y]
                SIG_ERR_BIN[Bin_i] = SIG_ERR_map[x][y]
                H3_ERR_BIN[Bin_i] = H3_ERR_map[x][y]
                H4_ERR_BIN[Bin_i] = H4_ERR_map[x][y]

LOS_vel = np.mean(VEL_BIN) #line of sight velocity

nump_data = np.transpose([VEL_BIN - LOS_vel,VEL_ERR_BIN,SIG_BIN,SIG_ERR_BIN,H3_BIN,H3_ERR_BIN,H4_BIN,H4_ERR_BIN])

np.savetxt('kinem_gh_SOME.txt', nump_data, fmt='%8.3f',
           header='v        v_err    sigma    sigma_err')
np.savetxt('bins_SOME.txt', BINS,fmt='%8.4f %8.4f %7d',
           header='Binning scheme file\nX_pix  Y_pix  binid')
