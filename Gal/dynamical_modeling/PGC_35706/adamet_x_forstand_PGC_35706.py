#!/usr/bin/python
'''
This program is an example of running observationally-constrained Schwarzschild models
(the FORSTAND code, Vasiliev&Valluri 2020).
It has several modes of operation:

1.  Run a model for a particular choice of parameters for the potential, orbit library, etc.
    (actually, a series of models with the same potential/orbit, but different mass-to-light ratios,
    in which the velocities are rescaled before fitting to observational constraints).
    Normally one would launch several copies of the script with different parameters (e.g. Mbh),
    the results will be collected in a summary text file,
    and each model's best-fit LOSVD and orbit library will be stored in separate files.

2.  Display an interactive plot with several panels:
  - kinematic maps (v,sigma,Gauss-Hermite moments) of the data or the model(s),
    one may click on any aperture and examine the LOSVD profile in both the data and the current model;
  - a 2d plot of chi2 values as a function of potential parameters for a grid of models,
    one may choose a model from the grid, load its LOSVD and display corresponding kinematic maps;
  - LOSVD in the given aperture (both data constraints with uncertainties and the current model);
  - distribution of orbit weights of the current model in two projections of integral space:
    mean radius vs. normalized squared total angular momentum [L/Lcirc(E)]^2, or
    mean radius vs. orbit inclination Lz/L; the size of symbols indicates orbit weights.

3.  Display some diagnostic plots useful for assessing the overall setup before running any models:
    projected and 3d density profiles along major and minor axes, location and masses of spatial bins,
    circular-velocity curve of the potential, and the observed values of v0 and sigma.
    The surface density and especially the deprojected 3d density help to check the photometric model;
    they should be smooth and have a sensible shape. In the 3d density plot, we also show the nodes of
    the spatial discretization grid, which should cover the region of interest, especially the central
    cusp (if we are interested in measuring Mbh, we need to have at least a few bins in the density
    profile within the expected radius of influence).
    The top right panel shows the values of density constraints (essentially, cell masses);
    ideally they should be quite uniform (with the exception of the innermost ones, which may be
    smaller if the grid is deliberately made finer in the central region). If the dynamical range
    spans more than 3-4 orders of magnitude, there is a risk that the smallest bins don't get any
    suitable orbits passing through them, so these density constraints would be impossible to satisfy,
    and the model either will be infeasible (if the density tolerance is set to zero) or will have
    a large and unpredictably varying penalty for violating these density constraints, which is also
    bad. In these cases one would need to adjust the grid parameters or even change the density
    discretization scheme to a different kind (e.g. classic instead of cylindrical, or vice versa).
    The bottom right panel shows the circular-velocity curve (split between mass components)
    for the model with the currently chosen parameters (Mbh, stellar M/L, etc.). For comparison,
    the values of kinematic constraints (v0 and sigma) in all bins are plotted against radius.
    In a thin and cold disk observed edge-on, the maximum value of |v0| should be close to the
    amplitude of the total circular velocity, but in general it will be smaller by some factor.
    This plot may be used to get a rough idea about the expected M/L: the amplitude of Vcirc in
    the potential is scaled by sqrt(Upsilon), where Upsilon (provided as a command-line argument)
    is the starting value for scanning the M/L axis automatically performed by the code.

4.  Prepare mock data for running the first two tasks.
    For this, one needs an N-body snapshot - in this example it should be produced by running
    a separate script  example_self_consistent_model3.py, which creates a three-component galaxy model
    (disk + bulge + DM halo) with a central black hole.
    It should be straightforward to feed in any other N-body snapshot, adjusting a few parameters.

This script uses various routines from the submodule agama.schwarzlib, which are general enough
to be used in any model fitting task. The specifics of each particular galaxy are encoded as numerous
parameters scattered througout this file.
Almost all needed parameters have reasonable default values in this example, but of course these are
not necessarily optimal for other purposes.
When adapting this example to your particular dataset, you will need to modify/adjust many parameters;
those which likely need to be changed are marked by [REQ],
other parameters which may be adjusted but have reasonable default values are marked by [OPT].
It is convenient to assign default values at the beginning of the script, and optionally change
some of them by providing command-line arguments in the form  name=value.

To run a complete example of constructing a grid of Schwarzschild models, do the following:

1.  Construct an N-body snapshot which will be used for the mock data, by running
    example_self_consistent_model3.py
    Among other things, it will produce the two files model_disk_final, model_bulge_final,
    which together contain the stellar particles of the model galaxy.

2.  Prepare the mock photometric and kinematic datasets, by running
    example_forstand.py do=mock
    It requires two Python modules 'mgefit' and 'vorbin'.
    This will produce an MGE file with surface density, and two kinematic datasets
    (low-resolution covering a large part of the galaxy, and high-res for the central region).
    Kinematic data (LOSVDs) in ~200 Voronoi bins are provided in two alternative forms:
  - histogrammed LOSVDs with values and their uncertainties in ~15 bins across the velocity axis;
  - Gauss-Hermite moments of LOSVDs (v, sigma, h3, h4, h5, h6 and their uncertainties).
    The models could be run on either of these alternative datasets, the choice is controlled
    by the command-line argument  hist=[y/n]  (the mock data are always generated for both cases).

3.  Examine the model setup by running
    example_forstand.py do=test
    (this is less relevant for the mock dataset, but could be quite helpful when working
    with real data, to check if the parameters are sensible, or to diagnose possible problems).

4.  Now run several series of models with different values of Mbh:
    example_forstand.py do=run Mbh=...
    (of course, one may also adjust many other parameters, including hist=[true/false])
    Each series of models has the same gravitational potential, scaled by several values of M/L;
    the LOSVD and the orbit properties of each model series are written into separate .npz files,
    and the summary of all models (grid of parameters and chi2 values) are stored in a single file
    results***.txt (*** is either GH or Hist).
    For each series of models with a given potential, the one with a M/L that gives the lowest chi2
    is converted into an N-body representation and written into a separate text file.
    The true value of Mbh is 1e8, and M/L is 1; it makes sense to explore at least a few series of
    models with Mbh ranging from 0 to ~3-5 times the true value.

5.  Finally, one may explore the grid of models for all values of Mbh and M/L by running
    example_forstand.py do=plot [hist=... and other params]

When adapting this script to a particular galaxy with existing observational datasets,
start from step 4 (to make sure that the observed kinematic maps look reasonable and
geometric parameters of the model, such as viewing angles and grids, are correct),
and then go back to step 3 (run several series of models) before going to step 4 again.

This script is mainly tailored to axisymmetric systems, although the Forstand code is applicable
in a more general context (e.g., to rotating triaxial barred galaxies).
The main limitation is the lack of suitable deprojection methods: in this example we use
the Multi-Gaussian expansion to fit the 2d surface density profile of the N-body snapshot
and then deproject it into an ellipsoidally stratified 3d density profile.
In the case of triaxial systems, especially with radially varying axis ratios, this procedure
is much less reliable and may even fail to produce a reasonable deprojection if the actual
3d shape is not well described by ellipsoids.
For triaxial systems, there are two rather than one angle specifying the orientation:
inclination (beta) and the angle alpha between the major axis and the line of nodes,
and in the rotating case, the pattern speed Omega is also a free parameter.

The model parameters and corresponding chi^2 values are stored in a single file
resultsGH.txt (for Gauss-Hermite parametrization of the LOSVD) or
resultsHist.txt (for histograms), and the kinematic maps and orbit distribution of each series
of models with the same potential and varying M/L are stored in a separate .npz archive.
In the interactive plotting regime, the likelihood surface is shown as a function of two
model parameters (in this example, Mbh and M/L), but one may choose another pair of parameters
by providing different columns of the results file as "aval", "bval" arguments for
agama.schwarzlib.runPlot(...); the values of remaining parameters should then be fixed and
specified as command-line arguments.
'''

import sys, numpy, agama

############### parse parameters from command-line arguments or assign default values #############

arglist = []
for arg in sys.argv[1:]:
    nameval = arg.split('=')
    if len(nameval)!=2:
        raise ValueError('Command-line arguments should be in the form  name=value')
    arglist.append([nameval[0].upper(), nameval[1]])
args = dict(arglist)

distance  = float(args.get('DISTANCE', 117490))# [REQ] assumed distance [kpc]
arcsec2kpc= distance * numpy.pi / 648000        # conversion factor (number of kiloparsecs in one arcsecond)
agama.setUnits(mass=1, length=arcsec2kpc, velocity=1)  # [OPT] units: mass = 1 Msun, length = 1", velocity = 1 km/s
Mbh       = float(args.get('MBH', 1e7))         # [REQ] mass of the central black hole  [Msun]
Omega     = float(args.get('OMEGA', 0))         # [REQ] pattern speed (relevant only for non-axisymmetric models) [km/s/length_unit]
halotype  =       args.get('HALOTYPE', 'nfw')   # [OPT] halo type: 'LOG' or 'NFW'
Upsilon   = float(args.get('UPSILON', 15.0))     # [OPT] initial value of mass-to-light ratio in the search
multstep  = float(args.get('MULTSTEP', 1.1))   # [OPT] multiplicative step for increasing/decreasing Upsilon during grid search
numOrbits = int  (args.get('NUMORBITS', 5000)) # [OPT] number of orbit in the model (size of orbit library)
intTime   = float(args.get('INTTIME', 100.0))   # [OPT] integration time in units of orbital period
regul     = float(args.get('REGUL', 1. ))       # [OPT] regularization parameter (larger => more uniform orbit weight distribution in models)
incl      = float(args.get('INCL', 42.0))       # [REQ] inclination angle (0 is face-on, 90 is edge-on) [degrees]
beta      = incl * numpy.pi/180                 # same in radians
alpha_deg = float(args.get('ALPHA', 0))       # [REQ] azimuthal angle of viewing direction in the model coordinates (relevant only for non-axisym)
alpha     = alpha_deg * numpy.pi/180            # same in radians
degree    = int  (args.get('DEGREE', 2))        # [OPT] degree of B-splines  (0 means histograms, 2 or 3 is preferred)
symmetry  = 'a'                                 # [OPT] symmetry of the model ('s'pherical, 'a'xisymmetric, 't'riaxial)
addnoise  = bool (args.get('ADDNOISE', True))   # [OPT] whether to add a realistic amount of noise in generating mock datacubes
seed      = int  (args.get('SEED', 0))          # [OPT] random seed (different values will create different realizations of the initial conditions for the orbit library)
nbody     = int  (args.get('NBODY', 100000))    # [OPT] number of particles for the N-body representation of the best-fit model
nbodyFormat = args.get('NBODYFORMAT', 'text')   # [OPT] format for storing N-body snapshots (text/nemo/gadget)
command   = args.get('DO', '').upper()          # [REQ] operation mode: 'RUN' - run a model, 'PLOT' - show the model grid and maps, 'TEST' - show diagnostic plots, 'MOCK' - create mock maps
usehist   = args.get('HIST', 'n')[0] in 'yYtT1' # [OPT] whether to use LOSVD histograms as input (default 'no' is to use GH moments)
variant   = 'Hist' if usehist else 'GH'         # suffix for disinguishing runs using histogramed LOSVDs or GH moments
fileResult= 'results%s.txt' % variant           # [OPT] filename for collecting summary results for the entire model grid
numpy.random.seed(32)                           # [OPT] make things repeatable when generating mock data (*not* related to the seed for the orbit library)
numpy.set_printoptions(precision=4, linewidth=9999, suppress=True)

# In this example, we use the Multi-Gaussian Expansion to parametrize
# the surface density profile and deproject it into the 3d density profile,
# but the code works with any other choice of 3d density model.
filenameMGE = 'mge_PGC35706_z_legacy.txt'    # [REQ] file with parameters of the MGE model for the surface density profile (if MGE is used)

### common parameters for kinematic datasets (though in principle they may also differ between them)
gridv  = numpy.linspace(-250, 250, 46)  # [REQ] the grid in model velocity space (will be multiplied by sqrt(Upsilon) when comparing to data)
velpsf = 0.0                  # [OPT] velocity-space PSF (usually not needed, as the spectroscopic fits produce deconvolved LOSVDs)
# [OPT] define the degree and velocity grid for the observed LOSVD provided as histograms or (less likely) higher-degree B-splines;
# the following two lines are needed [REQ] only if the input is provided in the form of binned LOSVDs (usehist=True),
# but we also use these parameters to generate mock LOSVD histograms if command=='MOCK'
hist_degree = 0               # [OPT] B-spline degree for the observed LOSVDs (0 means histogram)
hist_gridv  = numpy.linspace(-400, 400, 17)  # [OPT] velocity grid for the observed LOSVDs (boundaries of velocity bins, not centers!)

### parameters of the 1st kinematic dataset
gamma1 = 0.0 * numpy.pi/180  # [REQ] CW rotation angle of the image-plane X axis relative to the line of nodes (=major axis for axisym.systems)
psf1   = 1.0                  # [REQ] width of the Gaussian PSF ( may use more than one component: [ [width1, weight1], [width2, weight2] ] )
kinemParams1 = dict(          # parameters passed to the constructor of the Target class
    type     = 'LOSVD',
    symmetry = symmetry,      # symmetry properties of the potential
    alpha    = alpha,         # two angles defining the orientation of the model
    beta     = beta,          # w.r.t. image plane (same for all kinematic datasets)
    gamma    = gamma1,        # third angle is the rotation of the image plane, may be different for each dataset
    psf      = psf1,          # spatial PSF
    velpsf   = velpsf,        # velocity-space PSF
    degree   = degree,        # parameters for the internal datacube represented by B-splines:
    gridv    = gridv,         # usually will be identical for all datasets (except gridx,gridy which is determined by apertures)
)
filenameVorBin1 = 'bins_PGC35706.txt' # [REQ] Voronoi binning scheme for this dataset
filenameHist1   = 'kinem_hist_i%.0f_lr.txt'   % incl # [REQ] histogrammed representation of observed LOSVDs
filenameGH1     = 'kinem_gh_PGC35706.txt'# [REQ] Gauss-Hermite parametrization of observed LOSVDs (usually only one of these two
"""
### same for the 2nd kinematic dataset [OPT] - may have only one dataset, or as many as needed
gamma2 = 115.0 * numpy.pi/180
psf2   = 1.0                  # in this case it's a high-resolution IFU datacube
kinemParams2 = dict(
    type     = 'LOSVD',
    symmetry = symmetry,
    alpha    = alpha,
    beta     = beta,
    gamma    = gamma2,
    psf      = psf2,
    velpsf   = velpsf,
    degree   = degree,
    gridv    = gridv,
)
filenameVorBin2 = 'bins_LEDA2220522.txt'
filenameHist2   = 'kinem_hist_i%.0f_hr.txt'   % incl
filenameGH2     = 'kinem_gh_LEDA2220522.txt'
"""

def lnprob_fun(pars0):

    command = 'RUN'

    vhalo     = float(pars0[0])       # [OPT] asymptotic (LOG) or peak (NFW) circular velocity of the halo [km/s]
    rhalo     = float(pars0[1])       # [OPT] core (LOG) or scale (NFW) radius of the halo [lenth_unit]

### assemble the datasets (Targets and constraints)
    datasets = []

### 0: photometry => 3d density profile and its discretization scheme for a density Target
    
# read the input MGE file, skipping the first three lines as comments, deproject it and construct the Density object.
# Instead of MGE, one may use any other parametric density profile, e.g. one or more Sersic components with parameters
# determined by photometric fitting software such as Galfit
    try:
        mge = numpy.loadtxt(filenameMGE)   # [REQ] file with MGE parametrization of surface density profile
    except:
        print('%s not found; you need to generate the mock data first, as explained at the beginning of this file' % filenameMGE)
        exit()

    densityStars = agama.schwarzlib.makeDensityMGE(mge, distance, arcsec2kpc, beta)
#densityStars = agama.Density(agama.Density('dens_disk'), agama.Density('dens_bulge'))  # true profiles of this mock dataset

### parameters for the density dataset
# the choice of discretization scheme depends on the morphological type of the galaxy being modelled:
# for disky systems, DensityCylindrical[TopHat/Linear] is preferred, either in the axisymmetric regime
# (mmax=0), or more generally with mmax>0;
# for spheroidal systems, DensityClassic[TopHat/Linear] or DensitySphHarm may be more suitable,
# and in any case, the choice of radial [and possibly vertical] grid requires careful consideration.
# Here we do it automatically to ensure that the grid covers almost the entire model
# and has roughly equal mass in each shell (for spherical) or slab (for cylindrical grids),
# but this might not be suitable for every case; in particular, one may wish to make the grids
# denser in the central region to better constrain the 3d density profile near the black hole.
    densityParams = dict(type = (
        'DensityClassicTopHat',
        'DensityClassicLinear',
        'DensitySphHarm',
        'DensityCylindricalTopHat',
        'DensityCylindricalLinear')[2])   # [REQ] choose one of these types!
# use the discrete samples from the density profile to choose the grid parameters
    samples = densityStars.sample(10000)[0]
    if densityParams['type'] == 'DensityClassicTopHat' or densityParams['type'] == 'DensityClassicLinear':
    # create a grid in elliptical radius with axis ratio chosen to [roughly] match those of the density profile
        axes = numpy.percentile(numpy.abs(samples), 90, axis=0)  # three principal axes in the outer part of the profile
        axes/= numpy.exp(numpy.mean(numpy.log(axes)))  # normalize so that the product of three axes is unity
        ellrad = numpy.sum((samples / axes)**2, axis=1)**0.5
    # [OPT] make the inner grid segment contain 1% of the total mass
    # (to better constrain the density near the black hole, though this may need some further tuning),
    # and the remaining segments contain roughly equal fractions of mass up to 99% of the total mass
        densityParams['gridr'] = numpy.hstack([0, numpy.percentile(ellrad, tuple(numpy.linspace(1, 99, 24))) ])
        densityParams['axisRatioY'] = axes[1] / axes[0]
        densityParams['axisRatioZ'] = axes[2] / axes[0]
        print('%s grid in elliptical radius: %s, axis ratios: y/x=%.3g, z/x=%.3g' %
            (densityParams['type'], densityParams['gridr'], densityParams['axisRatioY'], densityParams['axisRatioZ']))
    # [OPT] each shell in the elliptical radius is divided in three equal 'panes'
    # adjacent to each of the principal axes, and then each pane is further divided
    # into a square grid of cells with stripsPerPane elements on each side
        densityParams['stripsPerPane'] = 2
    elif densityParams['type'] == 'DensitySphHarm':
        # this discretization scheme uses a grid in spherical radius and a spherical-harmonic expansion in angles
        sphrad = numpy.sum(samples**2, axis=1)**0.5
    # [OPT] same procedure as above, using roughly equal-mass bins in spherical radius except the innermost one
        densityParams['gridr'] = numpy.hstack([0, numpy.percentile(sphrad, tuple(numpy.linspace(1, 99, 24))) ])
        # [OPT] order of angular spherical-harmonic expansion in theta and phi (must be even)
        densityParams['lmax'] = 0 if symmetry[0]=='s' else 8
        densityParams['mmax'] = 0 if symmetry[0]!='t' else 6
        print('%s grid in spherical radius: %s, lmax=%i, mmax=%i' %
            (densityParams['type'], densityParams['gridr'], densityParams['lmax'], densityParams['mmax']))
    elif densityParams['type'] == 'DensityCylindricalTopHat' or densityParams['type'] == 'DensityCylindricalLinear':
        sampleR = (samples[:,0]**2 + samples[:,1]**2)**0.5
        samplez = abs(samples[:,2])
        # [OPT] choose the grids in R and z so that each 'slab' (1d projection along the complementary coordinate)
        # contains approximately equal mass, though this doesn't guarantee that the 2d cells would be even roughly balanced
        densityParams['gridR'] = numpy.hstack([0, numpy.percentile(sampleR, tuple(numpy.linspace(1, 99, 20))) ])
        densityParams['gridz'] = numpy.hstack([0, numpy.percentile(samplez, tuple(numpy.linspace(1, 99, 15))) ])
        # [OPT] number of azimuthal-harmonic coefficients (0 for axisymmetric systems)
        densityParams['mmax']  = 0 if symmetry[0]!='t' else 6
        print('%s grid in R: %s, z: %s, mmax=%i' %
            (densityParams['type'], densityParams['gridR'], densityParams['gridz'], densityParams['mmax']))

    datasets.append(agama.schwarzlib.DensityDataset(
        density=densityStars,
        # [OPT] fractional tolerance (e.g., 0.01) on the values of density constraints;
        # may be 0, requiring to satisfy them exactly, but in this case the solution may be infeasible
        tolerance=0.01,
        alpha=alpha,     # the orientation of intrinsic model coordinates w.r.t. the observed ones,
        beta=beta,       # specified by two Euler angles (used only for plotting the projected density)
        **densityParams  # remaining parameters set above
    ) )


### 1: 1st kinematic dataset

# read the Voronoi binning scheme and convert it to polygons (aperture boundaries)
    vorbin1    = numpy.loadtxt(filenameVorBin1)
    apertures1 = agama.schwarzlib.getBinnedApertures(xcoords=-vorbin1[:,0], ycoords=vorbin1[:,1], bintags=vorbin1[:,2])
# note that when using real observational data, the coordinate system in the image plane is usually
# right-handed, with Y pointing up and X pointing right. This is different from the convention used
# in Agama, where X points left. Therefore, one will need to invert the X axis of the observed dataset:
# getBinnedApertures(xcoords=-vorbin[:,0], ...)

# use either histograms or GH moments as input data
    if usehist:
        # [REQ] read the input kinematic data in the form of histograms;
        # if using the mock data as produced by this script, each line contains both values and errors for each velocity bin
        # in a given aperture, but when using data coming from other sources, may need to adjust the order of columns below
        kindat1 = numpy.loadtxt(filenameHist1)
        datasets.append(agama.schwarzlib.KinemDatasetHist(
            density   = densityStars,
            tolerance = 0.01,              # [REQ] relative error in fitting aperture mass constraints
            obs_val   = kindat1[:, 0::2],  # [REQ] values of velocity histograms
            obs_err   = kindat1[:, 1::2],  # [REQ] errors in these values
            obs_degree= hist_degree,
            obs_gridv = hist_gridv,
            apertures = apertures1,
            **kinemParams1
        ) )
    else:
        # [REQ] read the input kinematic data (V, sigma, higher Gauss-Hermite moments);
        # if using the mock data produced by this script, each line contains interleaved values and errors of v,sigma,h3...h6,
        # but when using data coming from other sources, may need to adjust the order of columns below
        kindat1 = numpy.loadtxt(filenameGH1)
        datasets.append(agama.schwarzlib.KinemDatasetGH(
            density   = densityStars,
            tolerance = 0.01,              # [REQ] relative error in fitting aperture mass constraints
            ghm_val   = kindat1[:, 0::2],  # [REQ] values of v,sigma,h3,h4...
            ghm_err   = kindat1[:, 1::2],  # [REQ] errors in the same order
            apertures = apertures1,
            **kinemParams1
        ) )
    

    # create a dark halo according to the provided parameters (type, scale radius and circular velocity)
    if rhalo>0 and vhalo>0:
        if   halotype.upper() == 'LOG':
            densityHalo = agama.schwarzlib.makeDensityLogHalo(rhalo, vhalo)
        elif halotype.upper() == 'NFW':
            densityHalo = agama.schwarzlib.makeDensityNFWHalo(rhalo, vhalo)
        else:
            raise ValueError('Invalid halo type')
    else:
        densityHalo = agama.Density(type='Plummer', mass=0, scaleRadius=1)  # no halo
    
# additional density component for constructing the initial conditions:
# create more orbits at small radii to better resolve the kinematics around the central black hole
    densityExtra = agama.Density(type='Dehnen', scaleradius=1)

# fiducialMbh: Mbh used to construct initial conditions (may differ from Mbh used to integrate orbits;
# the idea is to keep fiducialMbh fixed between runs with different Mbh, so that the initial conditions
# for the orbit library are the same, compensating one source of noise in chi2 due to randomness)
    fiducialMbh = densityStars.totalMass() * 0.01

    # potential of the galaxy, excluding the central BH
    pot_gal   = agama.Potential(type='Multipole',
        density=agama.Density(densityStars, densityHalo),  # all density components together
        lmax=32,  # lmax is set to a large value to accurately represent a disky density profile
        mmax=0 if symmetry[0]!='t' else 6, gridSizeR=40)  # mmax>0 only for triaxial systems
    # potential of the central BH
    pot_bh    = agama.Potential(type='Plummer', mass=Mbh, scaleRadius=1e-4)
    # same for the fiducial BH
    pot_bhfidu= agama.Potential(type='Plummer', mass=fiducialMbh, scaleRadius=1e-4)
    # total potential of the model (used to integrate the orbits)
    pot_total = agama.Potential(pot_gal, pot_bh)
    # total potential used to generate initial conditions only
    pot_fidu  = agama.Potential(pot_gal, pot_bhfidu)
    
    
    ### finally, decide what to do
    if command == 'RUN':
    
        # prepare initial conditions - use the same total potential independent of the actual Mbh
        # [OPT]: choose the sampling method: isotropic IC drawn from Eddington DF are created by
        #   density.sample(numorbits, potential)
        # while IC with preferential rotation (for disky models) are constructed from axisymmetric Jeans eqns by
        #   density.sample(numorbits, potential, beta={0-0.5}, kappa={1 or -1, depending on sign of rotation})
        # Here we add together two sets of IC - the majority of orbits sampled with axisymmetric Jeans eqns,
        # plus a small fraction additionally sampled from the central region to improve coverage.
        # Different values of the 'seed' parameter will create initial conditions with different number of orbits,
        # which effectively makes a completely new random sample, and then the number is truncated back to numOrbits.
        ic = numpy.vstack((
            densityStars.sample(int(numOrbits*0.85)+seed, potential=pot_fidu, beta=0.3, kappa=1)[0][seed:],
            densityExtra.sample(int(numOrbits*0.15), potential=pot_fidu)[0] ))
    
    
        # launch the orbit library and perform fits for several values of Upsilon;
        bestfit_chi2 = agama.schwarzlib.runModel(datasets=datasets, potential=pot_total, ic=ic,
            intTime=intTime, Upsilon=Upsilon, multstep=multstep, regul=regul, Omega=Omega,
            # [OPT] prefix - common part of the file name storing LOSVDs of each model in this series;
            # the value of Upsilon is appended to each filename;  here one may adjust the string format or the list of parameters to store
            filePrefix = 'M%.3g_O%.3g_Rh%.3g_Vh%.3g_i%.0f_a%.0f_N%d_R%.2f_%s_' % 
                (Mbh, Omega, rhalo, vhalo, incl, alpha_deg, numOrbits, regul, variant) + densityParams['type'],
            # [OPT] data stored at the beginning of each line (= a separate model with a given Upsilon) in the results/summary file;
            # usually should contains the same parameters as in filePrefix, but separated by tabs.
            # Keep track of the order of parameters - when reading the results file in the plotting part of this script, the order should be the same.
            # After the linePrefix, each line in the result file will contain the value of Upsilon, values of chi2 for each dataset,
            # regularization penalty, and the name of the file with LOSVD of that model.
            linePrefix = '\t'.join([ '%.3g' % Mbh, '%.3g' % Omega, '%.3g' % rhalo, '%.3g' % vhalo,
                '%.0f' % incl, '%.0f' % alpha_deg, '%d' % numOrbits, '%.2f' % regul ]),
            # [OPT] results/summary file
            fileResult = fileResult,
            # [OPT] parameters for the N-body snapshot representing the best-fit model
            nbody = nbody, nbodyFormat = nbodyFormat )
            
    else:
        exit('Nothing to do!')
    return -0.5*bestfit_chi2

from adamet.adamet import adamet

def bestfit_adamet():

    RHalo,VHalo=150,180
    nstep=40
    pars0 = numpy.array([RHalo,VHalo])    # Starting guess
    #fargs = (RHalo,VHalo)   # Parameters to pass to the lnprob function
    sigpar = [20, 20]     # Order of magnitude of the uncertainties
    bounds = numpy.array([[80, 100], [200, 200]])  # [[min(RHalo), min(VHalo)], [max(RHalo), max(VHalo)]]
    
    pars, lnprob = adamet(lnprob_fun, pars0, sigpar, bounds, nstep,
                          nprint=nstep/3, labels=['RHalo', 'VHalo'])
    bestfit = pars[numpy.argmax(lnprob)]
    perc = numpy.percentile(pars, [15.86, 84.14], axis=0)
    sig_bestfit = numpy.squeeze(numpy.diff(perc, axis=0)/2)   # half of 68% interval
    print(f"RHalo = {bestfit[0]:0.2f} +/- {sig_bestfit[0]:0.2f}")
    print(f"VHalo = {bestfit[1]:0.2f} +/- {sig_bestfit[1]:0.2f}")
    print(pars)
    
bestfit_adamet()
