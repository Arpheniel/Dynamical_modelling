import mgefit
from mgefit.find_galaxy import find_galaxy
from mgefit.mge_fit_1d import mge_fit_1d
from mgefit.sectors_photometry import sectors_photometry
from mgefit.mge_print_contours import mge_print_contours
from mgefit.mge_fit_sectors_twist import mge_fit_sectors_twist
from mgefit.sectors_photometry_twist import sectors_photometry_twist
from mgefit.mge_print_contours_twist import mge_print_contours_twist
from mgefit.mge_fit_sectors import mge_fit_sectors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from os import path
import vorbin
import sys
import argparse
from vorbin.voronoi_2d_binning import voronoi_2d_binning


def dist_circle(xc, yc, s):
    """
    Returns an array in which the value of each element is its distance from
    a specified center. Useful for masking inside a circular aperture.

    The (xc, yc) coordinates are the ones one can read on the figure axes
    e.g. when plotting the result of my find_galaxy() procedure.

    """
    x, y = np.ogrid[-yc:s[0] - yc, -xc:s[1] - xc]  # note yc before xc
    rad = np.sqrt(x ** 2 + y ** 2)

    return rad


def surf_br(surf, M_sun):
    return 10 ** (-(surf - M_sun - 21.57) / 2.5)


def magnitudo(flx):
    return 22.5 - 2.5 * np.log10(flx)


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gal", type=str, default="None",
                    help="enter the galaxy name using -g or --gal")  # name of a galaxy
parser.add_argument("-ph", "--phot", type=str, default="None",
                    help="enter the name of the photometry file using -f or --ftp")  # name of file with photometry
parser.add_argument("-sp", "--spec", type=str, default="None",
                    help="enter the spectroscopy filename using -s or --spec")  # name of file with spectroscopy
parser.add_argument("-b", "--band", type=str, default="None",
                    help="enter the name of band using  -b or --band")  # name of band
parser.add_argument("-psf", "--psf", type=str, default="None",
                    help="enter a value of psf in arcsec  -psf or --psf")  # value of psf
parser.add_argument("-sc", "--scale", type=str, default="None",
                    help="enter a value of number of arcsec in pixel  -sc or --scale")  # value of number of arcsec in pixel
parser.add_argument("-MSun", "--mag", type=str, default="None",
                    help="enter a value of absolute magnitude of Sun in chosen band -MSun or --mag")  # value of absolute magnitude of Sun
args = parser.parse_args()

GAL_NAME = args.gal  # name of a galaxy
PHOT = args.phot  # name of file with photometry
KINEM = args.spec  # name of file with spectroscopy
BAND = args.band # name of band
PSF = float(args.psf) # value of psf
SCALE = float(args.scale) # value of number of arcsec in pixel
M_Sun = float(args.mag) # value of absolute magnitude of Sun

print(f'name of a galaxy: {GAL_NAME}')
print(f'name of file with fotometry: {PHOT}')
print(f'name of file with spectroscopy: {KINEM}')
print(f'name of band: {BAND}')
print(f'psf value: {PSF}')
print(f'number of arcsec in pixel: {SCALE}')

file = fits.open(PHOT)  # open legacy image
ph_gal = file[1].data
path = PHOT.split('/')[0]

#scale = 0.262  # arcsec/pixel
#sigmapsf_g = 1.44676  # mean psf_g
#sigmapsf_r = 1.21627  # mean psf_r
#sigmapsf_z = 1.17582  # mean psf_g

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 6))
ax.imshow(ph_gal)
ax.set_title(BAND)
fig.savefig(path + '/' + f'{GAL_NAME}_{BAND}',format = "eps")
plt.clf()

fig = plt.figure(figsize=(9, 6))
gal = find_galaxy(ph_gal, binning= 4, fraction=0.6, level=None,
                nblob=1, plot=True, quiet=False)
fig.savefig(path + '/' + f'{GAL_NAME}_{BAND}_find',format = "eps")
plt.clf()
print(gal.pa)
gal.pa = 50
fig = plt.figure(figsize=(9, 6))
s_gal = sectors_photometry(ph_gal, gal.eps, gal.pa, gal.xpeak, gal.ypeak, minlevel=0, plot=1)
fig.savefig(path + '/' + f'{GAL_NAME}_{BAND}_secphot',format = "eps")

mge_gal = mge_fit_sectors(s_gal.radius, s_gal.angle, s_gal.counts, gal.eps,
                        bulge_disk=False, fignum=1, linear=False, negative=False,
                        ngauss=20, normpsf=1., outer_slope=4, plot=True, qbounds=None,
                        quiet=False, rbounds=None, scale=SCALE, sigmapsf=PSF, sol=0)

plt.savefig(path + '/' + f'{GAL_NAME}_{BAND}-1d_MGEfit',format = "eps")
total_counts, sigma, q_obs = mge_gal.sol  # assign the solution to variables
print(mge_gal.sol.T)  # Print a table of best-fitting MGE parameters
surf = total_counts / (2 * np.pi * q_obs * sigma ** 2) * SCALE ** 2  # the Gaussians peak surface brightness
print(f'Surface brightness of Gaussians peak, Flux_{BAND}/arcsec^2', surf)
surf = magnitudo(surf)
print(f'Surface brightness of Gaussians peak, Mag_{BAND}/arcsec^2', surf)
L_pc2 = surf_br(surf, M_sun=M_Sun)
print(f'Surface brightness of Gaussians peak, L_sun/pc^2, {BAND}', L_pc2)

nump_data = np.transpose(np.array([np.round(L_pc2, 3), np.round(sigma * SCALE, 3), np.round(q_obs, 3)]))
np.savetxt(path + '/' + f'mge_{GAL_NAME}_{BAND}_legacy.txt', nump_data, fmt='%12.6g %11.3f %11.4f',
            header='MGE file\nsurface_density  width  axis_ratio\n[Lsun/pc^2]   [arcsec]')

fig = plt.figure(figsize=(9, 6))
cnt,image,model = mge_print_contours(ph_gal, gal.pa, gal.xpeak, gal.ypeak, mge_gal.sol, scale=SCALE,
    binning=4, sigmapsf=PSF)

np.savetxt(path + '/' + f'converted_fits_{GAL_NAME}_{BAND}.txt', image)

np.savetxt(path + '/' + f'model_{GAL_NAME}_{BAND}.txt', model)


fig.savefig(path + '/' + f'{GAL_NAME}_{BAND}_contours',format = "eps")


spectro = fits.open(KINEM)  # open MANGA results of analysis

binind = spectro['BINID'].data

x1, y1 = np.shape(binind[0])
X_pix, Y_pix, Bins = [], [], []
v, v_err = [], []
sig, sig_err = [], []
# h3, h3_err = [], []
# h4, h4_err = [], []
# h5, h5_err = [], []
# h6, h6_err = [], []

for i in range(x1):
    for j in range(y1):
        X_pix.append(i - x1/2)
        Y_pix.append(j - y1/2)
        Bins.append(binind[0, i, j])
        v.append(spectro['STELLAR_VEL'].data[i, j])
        v_err.append(spectro['STELLAR_VEL_IVAR'].data[i, j])
        sig.append(spectro['STELLAR_SIGMA'].data[i, j])
        sig_err.append(spectro['STELLAR_SIGMA_IVAR'].data[i, j])
        # h3.append()
        # h3_err.append()
        # h4.append()
        # h4_err.append()
        # h5.append()
        # h5_err.append()
        # h6.append()
        # h6_err.append()

filter = ((np.array([v]) != 0))[0]

v = np.array(v)[filter]
v_err = 1/np.array(v_err)[filter]
sig = np.array(sig)[filter]
sig_err = 1/np.array(sig_err)[filter]

X_pix = np.array(X_pix)[filter]
Y_pix = np.array(Y_pix)[filter]
Bins = np.array(Bins)[filter]

_, idx_vel = np.unique(v, axis=0, return_index=True)
v = v[np.sort(idx_vel)]
v_err = v_err[np.sort(idx_vel)]
sig = sig[np.sort(idx_vel)]
sig_err = sig_err[np.sort(idx_vel)]

_, idx_bin = np.unique(Bins, axis=0, return_index=True)
Bins_reindex = Bins[np.sort(idx_bin)]

kinem_reindex_dict = {}
for i in range(len(v)):
    kinem_reindex_dict[Bins_reindex[i]] = (v[i], v_err[i], sig[i], sig_err[i])

kinem_reindex_dict_sort = {}
v_sort = []
v_err_sort = []
sig_err_sort = []
sig_sort = []

for i in sorted(kinem_reindex_dict):
    kinem_reindex_dict_sort[i] = kinem_reindex_dict[i]
v_reindex = list(kinem_reindex_dict_sort.values())

for i in v_reindex:
    v_sort.append(i[0])
    v_err_sort.append(i[1])
    sig_sort.append(i[2])
    sig_err_sort.append(i[3])

v_sort = np.array(v_sort)
v_err_sort = np.array(v_err_sort)
sig_sort = np.array(sig_sort)
sig_err_sort = np.array(sig_err_sort)

nump_data = np.transpose(np.array([v_sort, v_err_sort, sig_sort, sig_err_sort]))

np.savetxt(path + '/' + f'kinem_gh_{GAL_NAME}.txt', nump_data, fmt='%8.3f',
           header='v        v_err    sigma    sigma_err')

vor = np.column_stack((X_pix, Y_pix, Bins))
#vor = sorted(vor, key=lambda vor:vor[2])
np.savetxt(path + '/' + f'bins_{GAL_NAME}.txt', vor, fmt='%8.4f %8.4f %7d',
           header='Binning scheme file\nX_pix  Y_pix  binid')



"""
Usage example:
    python3 mge_kinem.py -ph PGC21856/mosaic-00073567-PGC21856-z-CCD3-image.fits -sp PGC21856/manga-8138-6102-MAPS-VOR10-MILESHC-MASTARSSP.fits -g PGC21856 -b z -psf 1.17582 -sc 0.262 -MSun 3.89

"""
