#!/usr/bin/python
'''
Clean grid-point runner for Schwarzschild models of PGC 35706.
Derived from adamet_x_forstand_PGC_35706.py by stripping the AdaMet wrapper
(lnprob_fun / bestfit_adamet) and exposing RHALO and VHALO as command-line
arguments. All other parameters match the AdaMet run so that new points
are consistent with the existing 27-point chain in resultsGH.txt.

Usage:
    python run_forstand_grid_PGC_35706.py RHALO=107 VHALO=154
    python run_forstand_grid_PGC_35706.py RHALO=107 VHALO=154 NUMORBITS=40000 SAVE_ORBITS=yes

Output:
    - One line per Upsilon value appended to resultsGH.txt
      (same format as AdaMet run, so the files merge naturally)
    - A single .npz file with LOSVD and orbit weights per (RHALO, VHALO) point
    - Best-fit N-body snapshot as .nbody file
    - On stdout: ">>> DONE: RHALO=... VHALO=... bestfit_chi2=..."  <- for log parsing

Note on reproducibility:
    numpy.random.seed(32) is set at startup, so initial conditions are
    deterministic within a single process. HOWEVER, re-running a point
    that was already computed in the AdaMet chain will give a slightly
    different chi2 (typically ~20 units at N=10000), because the AdaMet
    run used a different PRNG state when lnprob_fun was called at that
    point. This is the expected noise floor from the finite orbit library.
'''

import sys, numpy, agama

############### parse parameters from command-line ####################
arglist = []
for arg in sys.argv[1:]:
    nameval = arg.split('=')
    if len(nameval) != 2:
        raise ValueError('Command-line arguments should be in the form  name=value')
    arglist.append([nameval[0].upper(), nameval[1]])
args = dict(arglist)

# --- mandatory grid coordinates ---
if 'RHALO' not in args or 'VHALO' not in args:
    raise ValueError('RHALO and VHALO must both be specified. '
                     'Example: python run_forstand_grid_PGC_35706.py RHALO=107 VHALO=154')
rhalo = float(args['RHALO'])
vhalo = float(args['VHALO'])

# --- all other parameters match the adamet run ---
distance    = float(args.get('DISTANCE', 117490))
arcsec2kpc  = distance * numpy.pi / 648000
agama.setUnits(mass=1, length=arcsec2kpc, velocity=1)
Mbh         = float(args.get('MBH', 1e7))
Omega       = float(args.get('OMEGA', 0))
halotype    =       args.get('HALOTYPE', 'nfw')
Upsilon     = float(args.get('UPSILON', 15.0))
multstep    = float(args.get('MULTSTEP', 1.1))
numOrbits   = int  (args.get('NUMORBITS', 10000))
intTime     = float(args.get('INTTIME', 100.0))
regul       = float(args.get('REGUL', 0.0))
incl        = float(args.get('INCL', 42.0))
beta        = incl * numpy.pi/180
alpha_deg   = float(args.get('ALPHA', 0))
alpha       = alpha_deg * numpy.pi/180
degree      = int  (args.get('DEGREE', 2))
symmetry    = 'a'
seed        = int  (args.get('SEED', 0))
nbody       = int  (args.get('NBODY', 100000))
nbodyFormat =       args.get('NBODYFORMAT', 'text')
usehist     = args.get('HIST', 'n')[0] in 'yYtT1'
save_orb    = args.get('SAVE_ORBITS', 'n')[0] in 'yYtT1'
save_orbits_to = args.get('SAVE_ORBITS_TO', '.')

variant     = 'Hist' if usehist else 'GH'
fileResult  = args.get('FILERESULT', 'results%s.txt' % variant)
numpy.random.seed(32)
numpy.set_printoptions(precision=4, linewidth=9999, suppress=True)

############### static inputs (match adamet script) ###################
filenameMGE     = 'mge_PGC35706_z_legacy.txt'
filenameVorBin1 = 'bins_PGC35706_Damirs.txt'
filenameHist1   = 'kinem_hist_i%.0f_lr.txt' % incl
filenameGH1     = 'kinem_gh_PGC35706_Damirs.txt'

gridv       = numpy.linspace(-250, 250, 46)
velpsf      = 0.0
hist_degree = 0
hist_gridv  = numpy.linspace(-400, 400, 17)

gamma1 = 0.0 * numpy.pi/180
psf1   = 1.0
kinemParams1 = dict(
    type='LOSVD', symmetry=symmetry, alpha=alpha, beta=beta, gamma=gamma1,
    psf=psf1, velpsf=velpsf, degree=degree, gridv=gridv,
)

############### assemble datasets (photometry + kinematics) ############
datasets = []

try:
    mge = numpy.loadtxt(filenameMGE)
except:
    print('%s not found' % filenameMGE)
    sys.exit(1)

densityStars = agama.schwarzlib.makeDensityMGE(mge, distance, arcsec2kpc, beta)

# DensitySphHarm discretisation — same as in the adamet run
densityParams = dict(type='DensitySphHarm')
samples = densityStars.sample(10000)[0]
sphrad  = numpy.sum(samples**2, axis=1)**0.5
densityParams['gridr'] = numpy.hstack([0, numpy.percentile(sphrad, tuple(numpy.linspace(1, 99, 24)))])
densityParams['lmax']  = 0 if symmetry[0] == 's' else 8
densityParams['mmax']  = 0 if symmetry[0] != 't' else 6
print('%s grid in spherical radius: %s, lmax=%i, mmax=%i' %
      (densityParams['type'], densityParams['gridr'],
       densityParams['lmax'], densityParams['mmax']))

datasets.append(agama.schwarzlib.DensityDataset(
    density=densityStars, tolerance=0.01,
    alpha=alpha, beta=beta, **densityParams))

# kinematics
vorbin1    = numpy.loadtxt(filenameVorBin1)
apertures1 = agama.schwarzlib.getBinnedApertures(
    xcoords=-vorbin1[:,0], ycoords=vorbin1[:,1], bintags=vorbin1[:,2])

if usehist:
    kindat1 = numpy.loadtxt(filenameHist1)
    datasets.append(agama.schwarzlib.KinemDatasetHist(
        density=densityStars, tolerance=0.01,
        obs_val=kindat1[:, 0::2], obs_err=kindat1[:, 1::2],
        obs_degree=hist_degree, obs_gridv=hist_gridv,
        apertures=apertures1, **kinemParams1))
else:
    kindat1 = numpy.loadtxt(filenameGH1)
    datasets.append(agama.schwarzlib.KinemDatasetGH(
        density=densityStars, tolerance=0.01,
        ghm_val=kindat1[:, 0::2], ghm_err=kindat1[:, 1::2],
        apertures=apertures1, **kinemParams1))

############### gravitational potential ################################
if rhalo > 0 and vhalo > 0:
    if halotype.upper() == 'LOG':
        densityHalo = agama.schwarzlib.makeDensityLogHalo(rhalo, vhalo)
    elif halotype.upper() == 'NFW':
        densityHalo = agama.schwarzlib.makeDensityNFWHalo(rhalo, vhalo)
    else:
        raise ValueError('Invalid halo type')
else:
    densityHalo = agama.Density(type='Plummer', mass=0, scaleRadius=1)

densityExtra = agama.Density(type='Dehnen', scaleradius=1)
fiducialMbh  = densityStars.totalMass() * 0.01

pot_gal    = agama.Potential(type='Multipole',
    density=agama.Density(densityStars, densityHalo),
    lmax=32, mmax=0 if symmetry[0] != 't' else 6, gridSizeR=40)
pot_bh     = agama.Potential(type='Plummer', mass=Mbh,         scaleRadius=1e-4)
pot_bhfidu = agama.Potential(type='Plummer', mass=fiducialMbh, scaleRadius=1e-4)
pot_total  = agama.Potential(pot_gal, pot_bh)
pot_fidu   = agama.Potential(pot_gal, pot_bhfidu)

############### build orbit library initial conditions #################
ic = numpy.vstack((
    densityStars.sample(int(numOrbits*0.85) + seed, potential=pot_fidu, beta=0.3, kappa=1)[0][seed:],
    densityExtra.sample(int(numOrbits*0.15),         potential=pot_fidu)[0]))

############### run the Schwarzschild model ############################
print('>>> Running Schwarzschild model: '
      'RHALO=%.3g, VHALO=%.3g, MBH=%.3g, N=%d, regul=%.2f, save_orb=%s' %
      (rhalo, vhalo, Mbh, numOrbits, regul, save_orb))

bestfit_chi2 = agama.schwarzlib.runModel(
    datasets=datasets, potential=pot_total, ic=ic,
    intTime=intTime, Upsilon=Upsilon, multstep=multstep, regul=regul, Omega=Omega,
    filePrefix='M%.3g_O%.3g_Rh%.3g_Vh%.3g_i%.0f_a%.0f_N%d_R%.2f_%s_' %
        (Mbh, Omega, rhalo, vhalo, incl, alpha_deg, numOrbits, regul, variant)
        + densityParams['type'],
    linePrefix='\t'.join(['%.3g' % Mbh, '%.3g' % Omega,
                          '%.3g' % rhalo, '%.3g' % vhalo,
                          '%.0f' % incl, '%.0f' % alpha_deg,
                          '%d' % numOrbits, '%.2f' % regul]),
    fileResult=fileResult,
    nbody=nbody, nbodyFormat=nbodyFormat,
    save_orbits=save_orb, save_orbits_to=save_orbits_to)

# parse-friendly final line
print('>>> DONE: RHALO=%.1f VHALO=%.1f bestfit_chi2=%.4f' %
      (rhalo, vhalo, bestfit_chi2))
