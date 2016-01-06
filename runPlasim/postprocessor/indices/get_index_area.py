import os, sys
import numpy as np
from netCDF4 import Dataset

# Index calculation functions definition
def getAreaRange(srcVar, lat, lon, WEIGHTS, rng, win):
    """ Calculate the Mean Surface Temperature. """
    nt = srcVar.shape[0]
    index = np.empty((nt,))
    ilat = (lat >= win[1]) & (lat <= win[3])
    ilon = (lon >= win[0]) & (lon <= win[2])
    WEIGHTS = WEIGHTS[ilat][:, ilon]
    for k in np.arange(nt):
        srcVarWin = srcVar[k][ilat][:, ilon]
        irng = (srcVarWin >= rng[0]) & (srcVarWin <= rng[1])
        index[k] = WEIGHTS[irng].sum()
    return index

S = float(sys.argv[1])
restartState = sys.argv[2]
indexChoice = sys.argv[3]
firstYear = int(sys.argv[4])
lastYear  = int(sys.argv[5])
yearsPerFile = int(sys.argv[6])
daysPerYear = int(sys.argv[7])

srcDir = '../w2sb/'

# Ranges definition
Tf = 271.25 # Freezing temperature in Kelvin
rngTf20 = [Tf - 20, Tf]
rngTf10 = [Tf - 10, Tf]
rngTf = [0., Tf]

# Windows definitions [xll, yll, xur, yur]
winGlobal = [0., -90., 360., 90.]
winNHemi = [0., 0., 360., 90.]
winNPole = [0., 60., 360., 90.]
winNTrop = [0., 0., 360., 30.]
winEq = [0., -15., 360., 15.]

# Indices definitions
areaBelowTf20Glob = ('tsa', 'tsa', getAreaRange, 'areabelowtf20glob', rngTf20, winGlobal)
areaBelowTf20NHemi = ('tsa', 'tsa', getAreaRange, 'areabelowtf20nhemi', rngTf20, winNHemi)
areaBelowTf10NHemi = ('tsa', 'tsa', getAreaRange, 'areabelowtf10nhemi', rngTf10, winNHemi)
areaBelowTfNHemi = ('tsa', 'tsa', getAreaRange, 'areabelowtfnhemi', rngTf, winNHemi)

# Index choice
if indexChoice == 'areabelowtf20glob':
    indexDef = areaBelowTf20Glob
elif indexChoice == 'areabelowtf20nhemi':
    indexDef = areaBelowTf20NHemi
elif indexChoice == 'areabelowtf10nhemi':
    indexDef = areaBelowTf10NHemi
elif indexChoice == 'areabelowtfnhemi':
    indexDef = areaBelowTfNHemi
else:
    sys.exit('Unknown index choice: %s!' % indexChoice)

# Grid definition
# Spherical harmonics = Legendre, Fourier
nlat = 32
nlon = nlat*2
(lat, latWeights) = np.polynomial.legendre.leggauss(nlat)
lat = -np.arcsin(lat) * 180 / np.pi
lon = np.linspace(0., 360., nlon, endpoint=False)
latWeights /= 2
lonWeights = np.ones((nlon,)) / nlon
(LONWEIGHTS, LATWEIGHTS) = np.meshgrid(lonWeights, latWeights)
WEIGHTS = LONWEIGHTS * LATWEIGHTS

# Allocate
nt = (lastYear - firstYear + 1) * daysPerYear
index = np.empty((nt,))
firstYearOfFile = np.arange(firstYear, lastYear, yearsPerFile)
for k in np.arange(firstYearOfFile.shape[0]):
    # Read netCDF dataset
    ncFilePath = '%s/%s_%s_%d_%05d_%05d.nc' \
        % (srcDir, indexDef[0], restartState, S,
           firstYearOfFile[k],
           firstYearOfFile[k] + yearsPerFile - 1)
    print 'Reading %s...' % ncFilePath
    dset = Dataset(ncFilePath, 'r')
    srcVar = dset.variables[indexDef[1]][:]

    # Process variable
    indexFile = indexDef[2](srcVar, lat, lon, WEIGHTS, indexDef[4], indexDef[5])
    index[k*yearsPerFile*daysPerYear:(k+1)*yearsPerFile*daysPerYear] \
        = indexFile

    # Close netCDF dataset
    dset.close()

# Save index
dstDir = '%s_%d/' % (restartState, S)
os.system('mkdir %s 2> /dev/null' % dstDir)
np.savetxt('%s/%s_%s_%d_%05d_%05d.txt' \
               % (dstDir, indexDef[3], restartState, S,
                  firstYear, lastYear), index)


