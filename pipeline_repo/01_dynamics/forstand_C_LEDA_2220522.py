#!/usr/bin/python
'''
forstand_C_LEDA_2220522.py
=================================
APPROACH C = 1-D mass-normalization scan at a COSMOLOGICALLY FIXED NFW
scale radius.  Sibling of forstand_grid_R1_LEDA_2220522.py (the 2-D
(r_halo, v_halo) prior-anchored grid); the two coexist side-by-side
with disjoint output names (this version's outputs end with _C1).

RATIONALE
---------
At MaNGA coverage (~1 R_e) the data constrain essentially ONE mass
degree of freedom; the halo scale radius r_s is not measured (V&V 2020).
So instead of scanning a 2-D (r_s, v_max) grid we FIX r_s at the
cosmological abundance-matching expectation and scan only the overall
halo normalization (v_halo), with Upsilon marginalized per cell by the
sqrt(Y) trick.  The output is a 1-D chi^2(v_halo) curve; Delta_chi^2
gives a 1-D confidence interval, and the best cell is run with orbits
saved -- the orbit distribution is the science deliverable.

For LEDA this is especially well-motivated: DynPop already finds its
halo cosmologically anchored (r_s/R_e ~ 10.6, r_s stable cyl-vs-sph to
~1.8%, free fit == cosmological-concentration fit), so fixing r_s here
simply adopts a value the data are consistent with.

HOW r_s IS SET (per galaxy, in-script, pure numpy)
--------------------------------------------------
    physical M*  =  catalog logM* + 2*log10(1/h)        (h-correction)
    M*  --(Moster+13 z=0 SHMR, bisection-inverted)-->   M200
    M200 --(Dutton & Maccio 2014 z=0 c-M)-->            c200
    r200 from M200 & rho_crit ; r_s = r200/c200 (kpc)
    r_s[arcsec] = r_s[kpc] / arcsec2kpc   (uses THIS script's distance)
    V_max prior from (M200, c200) -> centres the v_halo scan.
Override the fixed radius with  rhalo_fix=<arcsec>  on the command line.

INHERITED UNCHANGED FROM R1
---------------------------
    regul=1.0 (V&V O(1)); global sqrt(Y) rescaling of stars+halo+BH
    (V&V 2020 §2.3); numOrbits=20000; ic generation in pot_fidu; MGE
    files; density type (DensityCylindricalLinear for LEDA);
    inclination/PA; distance; kinematic data; Voronoi binning; Mbh=1e7;
    subprocess isolation + 3-attempt retry resilience.

OUTPUT NAMES (all suffixed _C1, fully disjoint from _R1)
--------------------------------------------------------
    resultsGH_C1.txt           grid_chi2_GH_C1.npz
    resultsGH_C1_chunk<N>.txt  grid_chi2_GH_C1_chunk<N>.npz
    resultsGH_C1_verify.txt    verify_GH_C1.npz

USAGE
-----
    bash launch_C_LEDA.sh
    (4 chunks; each is one worker over a slice of the v_halo points)

PHASES
------
    phase=all       (default if chunk=-1):  scan + verify + final, single process
    phase=grid                            :  only the 1-D scan
    phase=merge                           :  read chunk files, merge into master
    phase=verify                          :  global-rescale verification of scan minimum
    phase=final                           :  final orbit-saving model at master scan min
    phase=postgrid                        :  merge + verify + final (after parallel scan)
    phase=onecell                         :  internal phase used by subprocess isolation
'''

import os, sys, glob, subprocess, tempfile, time, traceback, numpy, agama

############### parse parameters from command-line arguments or assign default values #############

arglist = []
for arg in sys.argv[1:]:
    nameval = arg.split('=')
    if len(nameval)!=2:
        raise ValueError('Command-line arguments should be in the form  name=value')
    arglist.append([nameval[0].upper(), nameval[1]])
args = dict(arglist)

distance  = float(args.get("DISTANCE", 108643))
arcsec2kpc= distance * numpy.pi / 648000
agama.setUnits(mass=1, length=arcsec2kpc, velocity=1)
Mbh       = float(args.get('MBH', 1e7))
Omega     = float(args.get('OMEGA', 0))
halotype  =       args.get('HALOTYPE', 'nfw')
Upsilon   = float(args.get('UPSILON', 7.0)) 
multstep  = float(args.get('MULTSTEP', 1.1))
numOrbits = int  (args.get('NUMORBITS', 20000))
intTime   = float(args.get('INTTIME', 100.0))
regul     = float(args.get('REGUL', 1.0 ))   # R1: default 1.0 (was 0.0)
incl      = float(args.get('INCL', 27.0))
beta      = incl * numpy.pi/180
alpha_deg = float(args.get('ALPHA', 0))
alpha     = alpha_deg * numpy.pi/180
degree    = int  (args.get('DEGREE', 2))
symmetry  = 'a'
addnoise  = bool (args.get('ADDNOISE', True))
seed      = int  (args.get('SEED', 0))
nbody     = int  (args.get('NBODY', 100000))
nbodyFormat = args.get('NBODYFORMAT', 'text')
usehist   = args.get('HIST', 'n')[0] in 'yYtT1'

# ====================================================================================
#   APPROACH C: cosmological NFW scale-radius prior  (pure numpy, z=0)
# ====================================================================================
# Catalog stellar mass (MaNGA DynPop Sersic, in units of h^-2 Msun); the
# h-correction below converts it to physical Msun.  LEDA 2220522: logM* = 9.741.
LOGMSTAR_CAT = float(args.get('LOGMSTAR_CAT', 9.741))

# Planck-like cosmology, consistent with the DynPop distance scale.
_H0   = 67.7                       # km/s/Mpc
_hh   = _H0 / 100.0
_Gkpc = 4.30091e-6                 # kpc (km/s)^2 / Msun
_H0_kpc  = _H0 / 1000.0            # (km/s)/kpc
_rho_crit = 3.0*_H0_kpc**2 / (8.0*numpy.pi*_Gkpc)   # Msun/kpc^3

def _moster_mstar_over_m200(logM200):
    # Moster+13 z=0 stellar-to-halo mass relation -> M*/M200
    logM1, N, b, g = 11.590, 0.0351, 1.376, 0.608
    x = 10.0**logM200 / 10.0**logM1
    return 2.0*N / (x**(-b) + x**(g))

def _invert_moster(logMstar):
    # bisection for logM200 such that M200 * (M*/M200)(M200) = M*
    target = 10.0**logMstar
    lo, hi = 9.0, 15.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if 10.0**mid * _moster_mstar_over_m200(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)

def _dutton_maccio_c200(logM200):
    a, b = 0.905, -0.101            # Dutton & Maccio 2014, z=0
    return 10.0**(a + b*numpy.log10(10.0**logM200 * _hh / 1e12))

def nfw_cosmo_prior(logMstar_cat, arcsec2kpc_):
    """Return (rs_arcsec, vmax_kms) from abundance matching + c-M relation."""
    logMstar = logMstar_cat + 2.0*numpy.log10(1.0/_hh)   # h-correction -> physical
    logM200  = _invert_moster(logMstar)
    M200     = 10.0**logM200
    c200     = _dutton_maccio_c200(logM200)
    r200     = (M200 / (4.0/3.0*numpy.pi*200.0*_rho_crit))**(1.0/3.0)   # kpc
    rs_kpc   = r200 / c200
    V200     = (_Gkpc*M200/r200)**0.5
    mu_c     = numpy.log(1.0+c200) - c200/(1.0+c200)
    Vmax     = 0.465 * V200 * (c200/mu_c)**0.5
    return rs_kpc/arcsec2kpc_, Vmax, dict(logM200=logM200, c200=c200,
                                          r200=r200, rs_kpc=rs_kpc, V200=V200)

_rs_prior_arcsec, _vmax_prior, _prior_info = nfw_cosmo_prior(LOGMSTAR_CAT, arcsec2kpc)
# r_s is FIXED at this cosmological value (override with rhalo_fix=<arcsec>).
RHALO_FIX = float(args.get('RHALO_FIX', _rs_prior_arcsec))
# 1-D v_halo scan brackets the V_max prior; widens both ways to trace the
# normalization degeneracy (low end ~ maximal-stellar, high end ~ heavy halo).
NVHALO        = int  (args.get('NVHALO', 12))
VHALO_FRAC_LO = float(args.get('VHALO_FRAC_LO', 0.40))
VHALO_FRAC_HI = float(args.get('VHALO_FRAC_HI', 1.80))

# Which v_halo the final orbit-saving run uses.  chi^2 is flat-bottomed /
# degenerate at ~1 R_e, so the formal grid minimum is NOT a meaningful
# normalization here -- LEDA is finalized at the cosmological prior instead.
#   FINAL_BY=prior    (default) -> v_halo = grid point nearest the V_max prior
#   FINAL_BY=chi2min            -> global chi^2 minimum on the master grid
#   FINAL_VHALO=<v>            -> force this exact v_halo (overrides FINAL_BY)
FINAL_BY    =       args.get('FINAL_BY', 'prior').lower()
FINAL_VHALO = float(args.get('FINAL_VHALO', 0.0))

print('=== Approach-C cosmological halo prior ===')
print('  logM*_cat   = %.3f  -> logM200 = %.3f  c200 = %.2f' %
      (LOGMSTAR_CAT, _prior_info['logM200'], _prior_info['c200']))
print('  r_s (fixed) = %.2f"  (%.2f kpc)' % (RHALO_FIX, RHALO_FIX*arcsec2kpc))
print('  V_max prior = %.1f km/s   scan x[%.2f, %.2f], %d pts' %
      (_vmax_prior, VHALO_FRAC_LO, VHALO_FRAC_HI, NVHALO))
print('==========================================')

# ---- chunked-execution control ----
phase     =       args.get('PHASE', 'all').lower()
chunk     = int  (args.get('CHUNK',   -1))
nchunks   = int  (args.get('NCHUNKS', 1 ))
variant   = 'Hist' if usehist else 'GH'

# Approach C: all output names tagged _C1, fully disjoint from the _R1 grid.
if chunk >= 0:
    if not (0 <= chunk < nchunks):
        raise ValueError('chunk=%d must satisfy 0 <= chunk < nchunks=%d' % (chunk, nchunks))
    is_chunk_mode = True
    fileResult    = 'results%s_C1_chunk%d.txt'        % (variant, chunk)
    fileResultVer = 'results%s_C1_verify_chunk%d.txt' % (variant, chunk)
    GRID_FILE     = 'grid_chi2_%s_C1_chunk%d.npz'     % (variant, chunk)
    if phase == 'all':
        phase = 'grid'
else:
    is_chunk_mode = False
    fileResult    = 'results%s_C1.txt'        % variant
    fileResultVer = 'results%s_C1_verify.txt' % variant
    GRID_FILE     = 'grid_chi2_%s_C1.npz'     % variant
MASTER_GRID = 'grid_chi2_%s_C1.npz' % variant   # C1 canonical merged file

print('=== Run configuration ===')
print('  variant     = %s' % variant)
print('  phase       = %s' % phase)
print('  chunk mode  = %s%s' % (is_chunk_mode, (' (%d/%d)' % (chunk, nchunks)) if is_chunk_mode else ''))
print('  GRID_FILE   = %s' % GRID_FILE)
print('  fileResult  = %s' % fileResult)
print('  numOrbits   = %d' % numOrbits)
print('=========================')

save_orb       = False
save_orbits_to = ""
numpy.random.seed(32)
numpy.set_printoptions(precision=4, linewidth=9999, suppress=True)

filenameMGE = 'mge_LEDA_2220522_z_legacy.txt'

gridv  = numpy.linspace(-250, 250, 46)
velpsf = 0.0
hist_degree = 0
hist_gridv  = numpy.linspace(-400, 400, 17)
gamma1 = 40.0 * numpy.pi/180
psf1   = 1.0
kinemParams1 = dict(type='LOSVD', symmetry=symmetry, alpha=alpha, beta=beta,
                    gamma=gamma1, psf=psf1, velpsf=velpsf, degree=degree, gridv=gridv)
filenameVorBin1 = 'bins_LEDA_2220522_Damirs.txt'
filenameHist1   = 'kinem_hist_i%.0f_lr.txt'   % incl
filenameGH1     = 'kinem_gh_LEDA_2220522_Damirs.txt'


# ====================================================================================
#                                 LNPROB / MODEL CALL
# ====================================================================================

def lnprob_fun(pars0, bake_upsilon=1.0):
    """bake_upsilon=1.0 -> standard sqrt(Y)-scan;
       bake_upsilon!=1  -> Upsilon baked into stellar density, single solve at Upsilon=1."""
    vhalo = float(pars0[0])
    rhalo = float(pars0[1])
    baked = (bake_upsilon != 1.0)
    tag   = 'baked Y=%.3f' % bake_upsilon if baked else 'trick scan'

    try:
        mge = numpy.loadtxt(filenameMGE)
    except:
        print('%s not found' % filenameMGE);  exit()
    # R1: when baking Upsilon, rescale ALL mass components by Y (V&V 2020 §2.3):
    #   stars   : Sigma_0 -> Sigma_0 * Y      (column 0 of MGE)
    #   halo    : vhalo   -> vhalo   * sqrt(Y)   (NFW: M_halo ~ vhalo^2, so M*Y -> v*sqrt(Y))
    #   SMBH    : Mbh     -> Mbh     * Y
    # The old `verify_*` script rescaled only stars and produced a physically
    # inconsistent control (Δχ² ~ +1000 against trick in the diag test, vs ~+20
    # after this fix).
    if baked:
        mge = mge.copy()
        mge[:, 0] *= bake_upsilon
        vhalo_use = vhalo * (bake_upsilon ** 0.5)
        Mbh_use   = Mbh   * bake_upsilon
    else:
        vhalo_use = vhalo
        Mbh_use   = Mbh

    densityStars = agama.schwarzlib.makeDensityMGE(mge, distance, arcsec2kpc, beta)

    densityParams = dict(type=(
        'DensityClassicTopHat','DensityClassicLinear','DensitySphHarm',
        'DensityCylindricalTopHat','DensityCylindricalLinear')[4])
    samples = densityStars.sample(10000)[0]
    if densityParams['type'] in ('DensityClassicTopHat','DensityClassicLinear'):
        axes = numpy.percentile(numpy.abs(samples), 90, axis=0)
        axes/= numpy.exp(numpy.mean(numpy.log(axes)))
        ellrad = numpy.sum((samples / axes)**2, axis=1)**0.5
        densityParams['gridr'] = numpy.hstack([0, numpy.percentile(ellrad, tuple(numpy.linspace(1, 99, 24)))])
        densityParams['axisRatioY'] = axes[1]/axes[0]
        densityParams['axisRatioZ'] = axes[2]/axes[0]
        densityParams['stripsPerPane'] = 2
    elif densityParams['type'] == 'DensitySphHarm':
        sphrad = numpy.sum(samples**2, axis=1)**0.5
        densityParams['gridr'] = numpy.hstack([0, numpy.percentile(sphrad, tuple(numpy.linspace(1, 99, 24)))])
        densityParams['lmax']  = 0 if symmetry[0]=='s' else 8
        densityParams['mmax']  = 0 if symmetry[0]!='t' else 6
    elif densityParams['type'] in ('DensityCylindricalTopHat','DensityCylindricalLinear'):
        sampleR = (samples[:,0]**2 + samples[:,1]**2)**0.5
        samplez = abs(samples[:,2])
        densityParams['gridR'] = numpy.hstack([0, numpy.percentile(sampleR, tuple(numpy.linspace(1, 99, 20)))])
        densityParams['gridz'] = numpy.hstack([0, numpy.percentile(samplez, tuple(numpy.linspace(1, 99, 15)))])
        densityParams['mmax']  = 0 if symmetry[0]!='t' else 6

    datasets = [agama.schwarzlib.DensityDataset(
        density=densityStars, tolerance=0.01,
        alpha=alpha, beta=beta, **densityParams)]

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

    # R1: use rescaled vhalo_use / Mbh_use (== input vhalo/Mbh when not baked)
    if rhalo > 0 and vhalo_use > 0:
        if   halotype.upper() == 'LOG': densityHalo = agama.schwarzlib.makeDensityLogHalo(rhalo, vhalo_use)
        elif halotype.upper() == 'NFW': densityHalo = agama.schwarzlib.makeDensityNFWHalo(rhalo, vhalo_use)
        else: raise ValueError('Invalid halo type')
    else:
        densityHalo = agama.Density(type='Plummer', mass=0, scaleRadius=1)
    densityExtra = agama.Density(type='Dehnen', scaleradius=1)
    fiducialMbh  = densityStars.totalMass() * 0.01

    pot_gal    = agama.Potential(type='Multipole',
        density=agama.Density(densityStars, densityHalo),
        lmax=32, mmax=0 if symmetry[0]!='t' else 6, gridSizeR=40)
    pot_bh     = agama.Potential(type='Plummer', mass=Mbh_use,    scaleRadius=1e-4)
    pot_bhfidu = agama.Potential(type='Plummer', mass=fiducialMbh, scaleRadius=1e-4)
    pot_total  = agama.Potential(pot_gal, pot_bh)
    pot_fidu   = agama.Potential(pot_gal, pot_bhfidu)

    ic = numpy.vstack((
        densityStars.sample(int(numOrbits*0.85)+seed, potential=pot_fidu, beta=0.3, kappa=1)[0][seed:],
        densityExtra.sample(int(numOrbits*0.15), potential=pot_fidu)[0]))

    if baked:
        ups_to_use, delta_chi2, mstep_to_use = 1.0, -1.0, 1.001
        result_file = fileResultVer
        # R1: filename records the RESCALED vhalo, Mbh as actually used.
        file_prefix = 'VER_Y%.3f_M%.3g_O%.3g_Rh%.3g_Vh%.3g_i%.0f_a%.0f_N%d_R%.2f_%s_' % (
            bake_upsilon, Mbh_use, Omega, rhalo, vhalo_use, incl, alpha_deg, numOrbits, regul, variant
        ) + densityParams['type']
        line_prefix = '\t'.join([
            '%.3g' % Mbh_use, '%.3g' % Omega, '%.3g' % rhalo, '%.3g' % vhalo_use,
            '%.0f' % incl, '%.0f' % alpha_deg, '%d' % numOrbits, '%.2f' % regul,
            'Y_baked=%.3f' % bake_upsilon, 'Vh0=%.3g' % vhalo, 'Mbh0=%.3g' % Mbh ])
    else:
        ups_to_use, delta_chi2, mstep_to_use = Upsilon, 100.0, multstep
        result_file = fileResult
        file_prefix = 'M%.3g_O%.3g_Rh%.3g_Vh%.3g_i%.0f_a%.0f_N%d_R%.2f_%s_' % (
            Mbh, Omega, rhalo, vhalo, incl, alpha_deg, numOrbits, regul, variant
        ) + densityParams['type']
        line_prefix = '\t'.join([
            '%.3g' % Mbh, '%.3g' % Omega, '%.3g' % rhalo, '%.3g' % vhalo,
            '%.0f' % incl, '%.0f' % alpha_deg, '%d' % numOrbits, '%.2f' % regul])

    print('  [lnprob_fun] %s' % tag, flush=True)
    bestfit_chi2 = agama.schwarzlib.runModel(
        datasets=datasets, potential=pot_total, ic=ic,
        intTime=intTime, Upsilon=ups_to_use, multstep=mstep_to_use,
        deltaChi2=delta_chi2, regul=regul, Omega=Omega,
        filePrefix=file_prefix, linePrefix=line_prefix,
        fileResult=result_file,
        nbody=nbody, nbodyFormat=nbodyFormat,
        save_orbits=save_orb, save_orbits_to=save_orbits_to)
    return -0.5 * bestfit_chi2


# ====================================================================================
#                                  HELPERS
# ====================================================================================

def parse_results_for(rhalo_val, vhalo_val, results_file):
    if not os.path.exists(results_file):
        return []
    matches = []
    rhalo_str = '%.3g' % rhalo_val
    vhalo_str = '%.3g' % vhalo_val
    with open(results_file) as f:
        for ln in f:
            parts = ln.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            if parts[2] != rhalo_str or parts[3] != vhalo_str:
                continue
            try:
                ups = float(parts[8])
            except ValueError:
                continue
            chi_cols = []
            for p in parts[9:]:
                try:
                    chi_cols.append(float(p))
                except ValueError:
                    break
            if not chi_cols: continue
            if len(chi_cols) >= 2:
                chi2_no_regul   = sum(chi_cols[:-1])
                chi2_with_regul = sum(chi_cols)
            else:
                chi2_no_regul = chi2_with_regul = chi_cols[0]
            matches.append((ups, chi2_no_regul, chi2_with_regul))
    return matches


def best_upsilon_for(rhalo_val, vhalo_val, results_file):
    m = parse_results_for(rhalo_val, vhalo_val, results_file)
    if not m: return None, None, None
    m.sort(key=lambda x: x[1])
    return m[0]


# ====================================================================================
#                                 GRID SEARCH
# ====================================================================================

# 18-point rhalo grid: extends the original geomspace(60, 400, 11) with 7 new
# Approach C: r_s FIXED at the cosmological prior (single rhalo row);
# the scan runs over v_halo only -> a 1-D chi^2(v_halo) curve.
GRID_RHALO = numpy.array([RHALO_FIX])                                          # 1 pt: fixed r_s [arcsec]
GRID_VHALO = numpy.geomspace(VHALO_FRAC_LO*_vmax_prior,
                             VHALO_FRAC_HI*_vmax_prior, NVHALO)                # NVHALO pts around V_max prior


def _load_grid(path, nR, nV):
    """Load chi2, ups arrays from an npz file if it exists.
    Exact match -> load as is.  Otherwise -> NaN.

    R1: auto-embedding of a sub-block from an older smaller grid is
    DISABLED. chi^2 values from regul=0 runs are not directly comparable
    to regul=1 fits (they differ by O(10) systematically), so inheriting
    them would corrupt the master. Only files written by THIS R1 pipeline
    are accepted."""
    chi2 = numpy.full((nR, nV), numpy.nan)
    ups  = numpy.full((nR, nV), numpy.nan)
    if not os.path.exists(path):
        return chi2, ups
    try:
        d = numpy.load(path)
    except Exception as e:
        print('  [load_grid] could not read %s (%s)' % (path, e))
        return chi2, ups

    file_rhalo = d['rhalo']
    file_vhalo = d['vhalo']

    # exact match -- the only case we accept in R1
    if (file_rhalo.shape == GRID_RHALO.shape and file_vhalo.shape == GRID_VHALO.shape and
        numpy.allclose(file_rhalo, GRID_RHALO) and numpy.allclose(file_vhalo, GRID_VHALO)):
        chi2 = d['chi2'].copy()
        if 'upsilon' in d.files:
            ups = d['upsilon'].copy()
        return chi2, ups

    print('  [load_grid] %s has different grid shape; ignoring (R1 disables auto-migration)' % path)
    return chi2, ups


# R1 + REBALANCE: if CELLS_FILE is given (an explicit list of (i j) lines),
# work only on those cells. Otherwise fall back to stripe ownership.
_cells_file = args.get('CELLS_FILE', '').strip()
_explicit_cells = None
if _cells_file:
    _explicit_cells = set()
    try:
        with open(_cells_file) as _fh:
            for _ln in _fh:
                _parts = _ln.split()
                if len(_parts) >= 2:
                    try:
                        _explicit_cells.add((int(_parts[0]), int(_parts[1])))
                    except ValueError:
                        pass
        print('  [CELLS_FILE] %s: %d cells assigned to this worker' %
              (_cells_file, len(_explicit_cells)), flush=True)
    except Exception as _e:
        print('  [CELLS_FILE] could not read %s (%s); falling back to stripe rule'
              % (_cells_file, _e), flush=True)
        _explicit_cells = None

def is_my_cell(i, j):
    if _explicit_cells is not None:
        return (i, j) in _explicit_cells
    if not is_chunk_mode:
        return True
    k = i * len(GRID_VHALO) + j
    return (k % nchunks) == chunk


def _run_cell_in_subprocess(rhalo_val, vhalo_val, max_attempts=3):
    """Compute one cell of the grid in an isolated child Python process,
    so a C-level crash inside agama/OpenMP/numpy does not kill the parent
    chunk.  Retries up to `max_attempts` times on subprocess failure.

    R1: max_attempts bumped 2 -> 3, sleep between retries 2s -> 30s,
    per-attempt status print made more explicit. Background: random
    python-process crashes have been observed on this server; spacing
    retries gives the OS / NFS / shared memory a chance to settle.

    Returns (chi2, upsilon) -- both NaN if all attempts crashed.
    """
    # Build subprocess argv: forward all original args except phase and
    # the cell-specific overrides.  This keeps chunk/nchunks/numOrbits/etc
    # consistent so the subprocess writes to the right fileResult etc.
    _SKIP_KEYS = {'PHASE', 'RHALO', 'VHALO', 'OUTPUT'}
    forwarded  = ['%s=%s' % (k, v) for (k, v) in arglist if k not in _SKIP_KEYS]

    env = dict(os.environ)
    env['PYTHONFAULTHANDLER'] = '1'   # catch SIGSEGV/SIGABRT/SIGFPE/SIGBUS

    chi2_val = numpy.nan
    ups_val  = numpy.nan
    # New in this version: separate the EACCES / "could not launch" case
    # from a real python crash. Launch-failures mean the interpreter or
    # script is currently unreachable (NFS glitch is the usual culprit);
    # the cure is to wait, not to retry-fast.
    EACCES_PAUSE_SEC = int(os.environ.get('EACCES_PAUSE_SEC', '120'))
    all_launch_failed = True   # becomes False the moment any subprocess actually starts

    for attempt in range(1, max_attempts + 1):
        fd, result_path = tempfile.mkstemp(
            suffix='.npz', prefix='cell_%d_' % os.getpid(), dir='/tmp')
        os.close(fd)
        try: os.remove(result_path)
        except OSError: pass

        cmd = [sys.executable, '-u', os.path.abspath(__file__),
               'phase=onecell',
               'rhalo=%.10g' % rhalo_val,
               'vhalo=%.10g' % vhalo_val,
               'output=%s' % result_path] + forwarded

        print('    [subprocess attempt %d/%d]' % (attempt, max_attempts),
              flush=True)
        launch_failed = False
        try:
            rc = subprocess.call(cmd, env=env)
        except Exception as e:
            print('    [subprocess attempt %d/%d] could not launch: %s' %
                  (attempt, max_attempts, e), flush=True)
            rc = -1
            launch_failed = True
        else:
            all_launch_failed = False   # exec succeeded at least

        if rc == 0 and os.path.exists(result_path):
            try:
                r = numpy.load(result_path)
                chi2_val = float(r['chi2'])
                ups_val  = float(r['upsilon'])
                os.remove(result_path)
            except Exception as e:
                print('    [subprocess attempt %d/%d] result unreadable: %s' %
                      (attempt, max_attempts, e), flush=True)
                try: os.remove(result_path)
                except OSError: pass
                continue
            # Whether finite or NaN, the subprocess RAN to completion.
            if numpy.isfinite(chi2_val):
                print('    -> chi2=%.3f  best Y=%.3f' % (chi2_val, ups_val),
                      flush=True)
            else:
                print('    -> NaN (Python exception inside cell; not retrying)',
                      flush=True)
            return chi2_val, ups_val, False
        else:
            print('    [subprocess attempt %d/%d] CRASHED  exit=%d  '
                  '(no result file)' % (attempt, max_attempts, rc),
                  flush=True)
            try: os.remove(result_path)
            except OSError: pass
            if attempt < max_attempts:
                if launch_failed:
                    # OS/NFS unreachable: long pause -- 30 s is useless here,
                    # the file system needs minutes to come back.
                    print('    [subprocess attempt %d/%d] launch failure looks transient (NFS?); '
                          'sleeping %d s before retry' %
                          (attempt, max_attempts, EACCES_PAUSE_SEC), flush=True)
                    time.sleep(EACCES_PAUSE_SEC)
                else:
                    # In-process python crash: short pause is fine.
                    time.sleep(30)

    print('    -> NaN after %d crashed attempts' % max_attempts, flush=True)
    return numpy.nan, numpy.nan, all_launch_failed


def bestfit_grid():
    nR, nV = len(GRID_RHALO), len(GRID_VHALO)
    # priority: own chunk file (most current) > master (pre-chunked run heritage)
    chi2, ups = _load_grid(GRID_FILE, nR, nV)
    if is_chunk_mode and os.path.exists(MASTER_GRID) and MASTER_GRID != GRID_FILE:
        chi2_m, ups_m = _load_grid(MASTER_GRID, nR, nV)
        mask = numpy.isfinite(chi2_m) & ~numpy.isfinite(chi2)
        chi2[mask] = chi2_m[mask]
        ups [mask] = ups_m [mask]
        print('  [bestfit_grid] inherited %d cells from master' % int(mask.sum()))
    print('  [bestfit_grid] %d cells loaded as already done before this run' %
          int(numpy.sum(numpy.isfinite(chi2))))

    my_cells = [(i,j) for i in range(nR) for j in range(nV) if is_my_cell(i,j)]
    total    = len(my_cells)
    print('  [bestfit_grid] this worker owns %d cells\n' % total)

    # Track consecutive cells where ALL retry attempts failed to even
    # launch the subprocess (NFS / interpreter unreachable). If too
    # many in a row, abort the chunk with a special exit code so the
    # respawn-launcher picks it up; the issue is system-wide and
    # spinning through more cells just turns them all into NaN.
    consec_launch_fail = 0
    MAX_CONSECUTIVE_LAUNCH_FAIL = int(os.environ.get('MAX_CONSECUTIVE_LAUNCH_FAIL', '5'))

    # Separately track the timestamp of the last SUCCESSFUL (finite-chi2)
    # write. mtime on the npz file gets bumped on every save including NaN
    # writes (e.g. when looping through EACCES-failed cells), so it is a
    # bad proxy for "is real progress happening?". We persist this stamp in
    # the npz itself and show it in show_status_R1.sh.
    last_progress_t = 0.0
    if os.path.exists(GRID_FILE):
        try:
            _d = numpy.load(GRID_FILE)
            if 'last_progress_t' in _d.files:
                last_progress_t = float(_d['last_progress_t'])
        except Exception:
            pass

    count = 0
    for (i, j) in my_cells:
        count += 1
        rhalo_val = float(GRID_RHALO[i])
        vhalo_val = float(GRID_VHALO[j])
        head = '[%d/%d (chunk %s)]  cell (%d,%d)  rhalo=%6.2f" (%5.1f kpc)  vhalo=%5.1f km/s' % (
            count, total, str(chunk) if is_chunk_mode else 'serial',
            i, j, rhalo_val, rhalo_val*arcsec2kpc, vhalo_val)
        if numpy.isfinite(chi2[i, j]):
            print('%s   skipped (chi2=%.2f, Y=%.2f)' % (head, chi2[i, j], ups[i, j]))
            continue
        # Run this cell in an isolated subprocess so a C-level crash in
        # agama/OpenMP/numpy does NOT kill the chunk.
        print('\n%s   running...' % head, flush=True)
        chi2_val, ups_val, all_launch_failed = _run_cell_in_subprocess(
            rhalo_val, vhalo_val, max_attempts=3)   # R1: 3 attempts
        chi2[i, j] = chi2_val
        ups [i, j] = ups_val
        if numpy.isfinite(chi2_val):
            last_progress_t = time.time()
        numpy.savez(GRID_FILE,
                    rhalo=GRID_RHALO, vhalo=GRID_VHALO,
                    chi2=chi2, upsilon=ups,
                    arcsec2kpc=arcsec2kpc, distance=distance,
                    incl=incl, Mbh=Mbh, numOrbits=numOrbits,
                    last_progress_t=last_progress_t)

        # Consecutive launch-failure tracker: ANY successful exec resets it.
        if all_launch_failed:
            consec_launch_fail += 1
            if consec_launch_fail >= MAX_CONSECUTIVE_LAUNCH_FAIL:
                print('\n  [bestfit_grid] %d consecutive cells with all-launch-failure: '
                      'looks like %s is unreachable system-wide. '
                      'Exiting (rc=2) so the respawn-launcher can retry later.'
                      % (consec_launch_fail, sys.executable),
                      flush=True)
                sys.exit(2)
        else:
            consec_launch_fail = 0

    if numpy.all(numpy.isnan(chi2)):
        # R1: with REBALANCE, a chunk may legitimately end up with zero
        # assigned cells (all its work went elsewhere). Treat that as
        # clean completion, not a hard failure.
        if total == 0:
            print('  [bestfit_grid] no cells assigned to this worker; nothing to do.', flush=True)
            return numpy.array([float('nan'), float('nan')]), (0, 0)
        raise RuntimeError('No finite chi^2 values - this worker did nothing.')

    i_b, j_b   = numpy.unravel_index(numpy.nanargmin(chi2), chi2.shape)
    best_rhalo = float(GRID_RHALO[i_b])
    best_vhalo = float(GRID_VHALO[j_b])
    best_chi2  = float(chi2[i_b, j_b])
    best_ups   = float(ups [i_b, j_b]) if numpy.isfinite(ups[i_b, j_b]) else numpy.nan

    print('\n' + '='*70)
    print('Best known cell after this worker (over %d done):' %
          int(numpy.sum(numpy.isfinite(chi2))))
    print('  rhalo=%.2f"  vhalo=%.1f km/s  Y=%.3f  chi2=%.3f' %
          (best_rhalo, best_vhalo, best_ups, best_chi2))
    print('='*70)

    return numpy.array([best_vhalo, best_rhalo]), (i_b, j_b)


# ====================================================================================
#                                 MERGE
# ====================================================================================

def merge_chunks(nchunks):
    nR, nV = len(GRID_RHALO), len(GRID_VHALO)
    chi2, ups = _load_grid(MASTER_GRID, nR, nV)
    n_pre = int(numpy.sum(numpy.isfinite(chi2)))
    print('Master grid before merge: %d/%d cells' % (n_pre, nR*nV))

    chunk_files = sorted(glob.glob('grid_chi2_%s_C1_chunk*.npz' % variant))
    if not chunk_files:
        print('  No chunk files found.  Nothing to merge.')
    for f in chunk_files:
        ci, ui = _load_grid(f, nR, nV)
        mask = numpy.isfinite(ci) & ~numpy.isfinite(chi2)   # only fill NaN cells
        chi2[mask] = ci[mask]
        ups [mask] = ui[mask]
        print('  %s: added %d cells' % (f, int(mask.sum())))

    n_post = int(numpy.sum(numpy.isfinite(chi2)))
    numpy.savez(MASTER_GRID,
                rhalo=GRID_RHALO, vhalo=GRID_VHALO,
                chi2=chi2, upsilon=ups,
                arcsec2kpc=arcsec2kpc, distance=distance,
                incl=incl, Mbh=Mbh, numOrbits=numOrbits)
    print('Master grid after merge:  %d/%d cells   (gained %d)' % (n_post, nR*nV, n_post-n_pre))
    return chi2, ups


# ====================================================================================
#                              UPSILON-TRICK VERIFICATION
# ====================================================================================

def verify_upsilon_trick(i_best, j_best, max_cells=9):
    d = numpy.load(MASTER_GRID)
    rhalo_grid = d['rhalo']
    vhalo_grid = d['vhalo']
    chi2_grid  = d['chi2']
    ups_grid   = d['upsilon'] if 'upsilon' in d.files else numpy.full_like(chi2_grid, numpy.nan)

    candidates = []
    for di in range(-1, 2):
        for dj in range(-1, 2):
            i, j = i_best+di, j_best+dj
            if 0 <= i < chi2_grid.shape[0] and 0 <= j < chi2_grid.shape[1] and numpy.isfinite(chi2_grid[i, j]):
                candidates.append((i, j, chi2_grid[i, j]))
    candidates.sort(key=lambda t: t[2])
    candidates = candidates[:max_cells]

    print('\n' + '='*70)
    print('VERIFICATION of the sqrt(Y) trick: %d cells around the grid minimum' % len(candidates))
    print('='*70)

    results = []
    for k, (i, j, _) in enumerate(candidates):
        rh = float(rhalo_grid[i])
        vh = float(vhalo_grid[j])
        chi2_trick = float(chi2_grid[i, j])
        ups = float(ups_grid[i, j]) if numpy.isfinite(ups_grid[i, j]) else None
        if ups is None:
            # fall back to results files (master and chunks)
            for rf in [fileResult] + sorted(glob.glob('results%s_C1_chunk*.txt' % variant)):
                ups_tmp, _, _ = best_upsilon_for(rh, vh, rf)
                if ups_tmp is not None:
                    ups = ups_tmp; break
        marker = '*** MIN' if (i==i_best and j==j_best) else '       '
        print('\n%s [%d/%d] cell (%d,%d): rhalo=%6.2f"  vhalo=%5.1f km/s  Y=%s' %
              (marker, k+1, len(candidates), i, j, rh, vh,
               '%.3f' % ups if ups is not None else 'NaN'))
        print('              chi2(trick) = %.3f' % chi2_trick)
        if ups is None:
            print('              skipping: no Upsilon found anywhere'); continue
        try:
            lnp_h = lnprob_fun([vh, rh], bake_upsilon=ups)
            chi2_honest = -2.0 * lnp_h
            delta = chi2_honest - chi2_trick
            verdict = ('OK' if abs(delta)<1 else 'MILD' if abs(delta)<5 else
                       'BIASED' if abs(delta)<20 else 'STRONG-BIAS')
            print('              chi2(honest) = %.3f' % chi2_honest)
            print('              dchi2 = %+.3f  [%s]' % (delta, verdict))
        except Exception as e:
            chi2_honest = numpy.nan
            print('              FAILED: %s' % e)
        results.append((rh, vh, ups, chi2_trick, chi2_honest, i, j))

    if results:
        arr = numpy.array([(r[0],r[1],r[2],r[3],r[4]) for r in results],
                          dtype=[('rhalo','f8'),('vhalo','f8'),('upsilon','f8'),
                                 ('chi2_trick','f8'),('chi2_honest','f8')])
        numpy.savez('verify_%s_C1.npz' % variant, table=arr, i_best=i_best, j_best=j_best)
        print('\n%-8s %-8s %-7s %-12s %-12s %-10s' % (
            'rhalo"','vhalo','Y','chi2_trick','chi2_honest','dchi2'))
        print('-'*70)
        for r in results:
            d_ = r[4]-r[3] if numpy.isfinite(r[4]) else float('nan')
            print('%-8.2f %-8.1f %-7.3f %-12.3f %-12.3f %+10.3f' % (
                r[0], r[1], r[2], r[3], r[4], d_))
        deltas = numpy.array([r[4]-r[3] for r in results if numpy.isfinite(r[4])])
        if len(deltas):
            print('\nMedian |dchi2|: %.2f   Max |dchi2|: %.2f' %
                  (numpy.median(numpy.abs(deltas)), numpy.max(numpy.abs(deltas))))
    return results


# ====================================================================================
#                                       MAIN
# ====================================================================================

def assert_grid_complete():
    """Refuse to finalize on a grid with holes.  nanargmin silently skips
    NaN cells, so a crashed/Python-swallowed cell would let postgrid pick a
    wrong point without any warning.  This makes the missing cell fatal and
    prints exactly which v_halo points are absent."""
    if not os.path.exists(MASTER_GRID):
        raise RuntimeError('master grid %s not found; run phase=merge first' % MASTER_GRID)
    chi2 = numpy.load(MASTER_GRID)['chi2']
    bad  = [float(GRID_VHALO[j]) for i in range(chi2.shape[0])
                                 for j in range(chi2.shape[1])
                                 if not numpy.isfinite(chi2[i, j])]
    if bad:
        raise RuntimeError(
            'master grid incomplete: %d/%d cells missing (v_halo=%s). '
            'Re-run the scan (it skips done cells) before finalizing.'
            % (len(bad), chi2.size, ', '.join('%.1f' % v for v in bad)))
    print('  [assert_grid_complete] %d/%d cells finite -- OK' % (chi2.size, chi2.size))


def final_cell_from_master():
    """Pick the v_halo for the final orbit-saving run.
    Default = cosmological prior (chi^2 is degenerate here, see header)."""
    assert_grid_complete()
    d = numpy.load(MASTER_GRID)
    chi2 = d['chi2']
    if FINAL_VHALO > 0:
        j_b = int(numpy.argmin(numpy.abs(GRID_VHALO - FINAL_VHALO)))
        how = 'forced FINAL_VHALO=%.1f' % FINAL_VHALO
    elif FINAL_BY == 'chi2min':
        _, (i_b, j_b) = best_cell_from_master()
        how = 'chi2 minimum'
    else:  # 'prior'
        j_b = int(numpy.argmin(numpy.abs(GRID_VHALO - _vmax_prior)))
        how = 'cosmological prior V_max=%.1f' % _vmax_prior
    i_b = 0   # single fixed-r_s row
    bestfit = numpy.array([float(GRID_VHALO[j_b]), float(GRID_RHALO[i_b])])
    print('Final v_halo chosen by %s: (%d,%d)  rhalo=%.2f"  vhalo=%.1f km/s  chi2=%.3f' %
          (how, i_b, j_b, GRID_RHALO[i_b], GRID_VHALO[j_b], chi2[i_b, j_b]))
    return bestfit, (i_b, j_b)


def best_cell_from_master():
    if not os.path.exists(MASTER_GRID):
        raise RuntimeError('master grid %s not found; run phase=merge first' % MASTER_GRID)
    d = numpy.load(MASTER_GRID)
    chi2 = d['chi2']
    if numpy.all(numpy.isnan(chi2)):
        raise RuntimeError('master grid has no finite chi^2 cells')
    i_b, j_b = numpy.unravel_index(numpy.nanargmin(chi2), chi2.shape)
    bestfit  = numpy.array([float(GRID_VHALO[j_b]), float(GRID_RHALO[i_b])])
    print('Best cell on master grid: (%d,%d)  rhalo=%.2f"  vhalo=%.1f km/s  chi2=%.3f' %
          (i_b, j_b, GRID_RHALO[i_b], GRID_VHALO[j_b], chi2[i_b, j_b]))
    return bestfit, (i_b, j_b)


def main():
    global save_orb, save_orbits_to
    print('\n' + '#'*70)
    print('# PHASE: %s' % phase.upper())
    print('#'*70)

    if phase == 'all':
        bestfit, ij = bestfit_grid()
        # serial run: own GRID_FILE IS the master, so copy it for verify to find
        if not os.path.exists(MASTER_GRID):
            import shutil; shutil.copy(GRID_FILE, MASTER_GRID)
        verify_upsilon_trick(*ij)
        save_orb       = True
        save_orbits_to = "/data1/vgorad/dynam_mod/Dynamical_modelling/Gal/chemo-dynamical_modeling/LEDA_2220522"
        lnprob_fun(bestfit)
        print('\nDone.')

    elif phase == 'grid':
        bestfit, ij = bestfit_grid()
        print('\nWorker done.  Grid scan complete for this chunk.')

    elif phase == 'onecell':
        # Compute one cell, write {chi2, upsilon} to args['OUTPUT'], exit 0.
        # This phase is invoked by the parent (phase=grid) via subprocess.
        # A C-level crash inside lnprob_fun -> agama runModel will kill THIS
        # subprocess; the parent sees non-zero exit and retries.
        rhalo_one = float(args['RHALO'])
        vhalo_one = float(args['VHALO'])
        output    = args['OUTPUT']
        chi2_one  = numpy.nan
        ups_one   = numpy.nan
        try:
            lnp = lnprob_fun([vhalo_one, rhalo_one])
            chi2_one = -2.0 * lnp
            u_best, _, _ = best_upsilon_for(rhalo_one, vhalo_one, fileResult)
            ups_one  = float(u_best) if u_best is not None else numpy.nan
            print('  [onecell] chi2=%.3f  Y=%.3f' % (chi2_one, ups_one), flush=True)
        except Exception as e:
            # Python-level exception inside the model. Previously only str(e) was
            # printed, which hid WHERE it came from (no file/line). Print the FULL
            # traceback so the failing call in schwarzlib/agama is visible in the
            # chunk log, and drop a .failed marker next to the output so the parent
            # can tell a genuine crash apart from an honest NaN.
            print('  [onecell] Python exception: %r' % e, flush=True)
            traceback.print_exc()
            sys.stdout.flush(); sys.stderr.flush()
            try:
                with open(output + '.failed', 'w') as _fh:
                    _fh.write(traceback.format_exc())
            except Exception:
                pass
        numpy.savez(output, chi2=chi2_one, upsilon=ups_one)
        sys.exit(0)

    elif phase == 'merge':
        merge_chunks(nchunks)

    elif phase == 'verify':
        if not os.path.exists(MASTER_GRID) or not numpy.any(numpy.isfinite(numpy.load(MASTER_GRID)['chi2'])):
            merge_chunks(nchunks)
        bestfit, ij = best_cell_from_master()
        verify_upsilon_trick(*ij)

    elif phase == 'final':
        bestfit, ij = final_cell_from_master()
        save_orb       = True
        save_orbits_to = "/data1/vgorad/dynam_mod/Dynamical_modelling/Gal/chemo-dynamical_modeling/LEDA_2220522"
        lnprob_fun(bestfit)
        print('\nFinal orbit-saving run complete.')

    elif phase == 'postgrid':
        merge_chunks(nchunks)
        assert_grid_complete()                 # refuse to finalize on a holey grid
        _, ij_min = best_cell_from_master()    # verify the sqrt(Y) trick at the chi2 min
        verify_upsilon_trick(*ij_min)
        bestfit, ij = final_cell_from_master() # but finalize at the prior (degenerate chi2)
        save_orb       = True
        save_orbits_to = "/data1/vgorad/dynam_mod/Dynamical_modelling/Gal/chemo-dynamical_modeling/LEDA_2220522"
        lnprob_fun(bestfit)
        print('\nPost-grid pipeline complete.')

    else:
        raise ValueError('Unknown phase: %s' % phase)


if __name__ == '__main__':
    main()
