import os, sys
import numpy as np
from netCDF4 import Dataset

# Mean radius or Earth in meters
Rm = 6371000

# Index calculation functions definition
def getMean(srcVar, lat, lon, WEIGHTS, win1, win2):
    """ Calculate the Mean Surface Temperature. """
    nt = srcVar.shape[0]
    index = np.empty((nt,))
    ilat = (lat >= win1[1]) & (lat <= win1[3])
    ilon = (lon >= win1[0]) & (lon <= win1[2])
    WEIGHTS = WEIGHTS[ilat][:, ilon]
    WEIGHTS /= WEIGHTS.sum()
    for k in np.arange(nt):
        index[k] = (srcVar[k][ilat][:, ilon] * WEIGHTS).sum()
    return index

def getMerGradient(srcVar, lat, lon, WEIGHTS, win1, win2):
    """ Calculate the Zonal Mean Meridional Temperature Gradient """
    nt = srcVar.shape[0]
    index = np.empty((nt,))
    ilat1 = (lat >= win1[1]) & (lat <= win1[3])
    ilon1 = (lon >= win1[0]) & (lon <= win1[2])
    ilat2 = (lat >= win2[1]) & (lat <= win2[3])
    ilon2 = (lon >= win2[0]) & (lon <= win2[2])
    WEIGHTS1 = WEIGHTS[ilat1][:, ilon1]
    WEIGHTS1 /= WEIGHTS1.sum()
    WEIGHTS2 = WEIGHTS[ilat2][:, ilon2]
    WEIGHTS2 /= WEIGHTS2.sum()
    for k in np.arange(nt):
        index[k] = (srcVar[k][ilat2][:, ilon2] * WEIGHTS2).sum() \
            - (srcVar[k][ilat1][:, ilon1] * WEIGHTS1).sum()
    return index

def getEnergyTransport(srcVar, lat, lon, WEIGHTS, win1, win2):
    """ Calculate the Mean Surface Temperature. """
    nt = srcVar.shape[0]
    index = np.empty((nt,))
    for k in np.arange(nt):
        index[k] = (np.abs(srcVar[k]) * WEIGHTS).sum() / 4
    return index

def getLatMaxTransport(srcVar, lat, lon, WEIGHTS, win1, win2):
    """ Calculate the Mean Surface Temperature. """
    nt = srcVar.shape[0]
    index = np.empty((nt,))
    ilat = (lat >= win1[1]) & (lat <= win1[3])
    ilon = (lon >= win1[0]) & (lon <= win1[2])
    WEIGHTS = WEIGHTS[ilat][:, ilon]
    WEIGHTS /= WEIGHTS.sum()
    lat = lat[ilat]
    for k in np.arange(nt):
        # Get zonal averages of the NTR for the atmosphere
        indexLat = (srcVar[k][ilat][:, ilon] * WEIGHTS).sum(1)
        # Get the gradient
        gradLag = -(indexLat[:-1] - indexLat[1:]) / (lat[:-1] - lat[1:]) / Rm
        # Take the maximum
        index[k] = lat[np.argmax(gradLag)]
    return index


S = float(sys.argv[1])
restartState = sys.argv[2]
indexChoice = sys.argv[3]
firstYear = int(sys.argv[4])
lastYear  = int(sys.argv[5])
yearsPerFile = int(sys.argv[6])
daysPerYear = int(sys.argv[7])

# S = 1260
# restartState = 'cold'
# firstYear = 101
# lastYear = 1100
# yearsPerFile = 100
# daysPerYear = 360
# indexChoice = 'globmst'


srcDir = '../w2sb/'
# Windows definitions [xll, yll, xur, yur]
winGlob = [0., -90., 360., 90.]
winNHemi = [0., 0., 360., 90.]
winNPole = [0., 60., 360., 90.]
winNTrop = [0., 0., 360., 30.]
winEq = [0., -15., 360., 15.]
winEqL = [0., -30., 360., 30.]

# Indices definitions
GlobMST = ('tsa', 'tsa', getMean, 'globmst', winGlob, winGlob) # Mean Surface Temperature
NPoleMST = ('tsa', 'tsa', getMean, 'npolemst', winNPole, winNPole)
NTropMST = ('tsa', 'tsa', getMean, 'ntropmst', winNTrop, winNTrop)
EqMST = ('tsa', 'tsa', getMean, 'eqmst', winEq, winEq)

NPoleNTR = ('ntr', 'ntr', getMean, 'npolentr', winNPole, winNPole) # Mean Net TOA Radiative fluxes
EqLNTR = ('ntr', 'ntr', getMean, 'eqlntr', winEqL, winEqL) # Mean Net TOA Radiative fluxes
EnergyTransport = ('ntr', 'ntr', getEnergyTransport, 'energytransport', winEqL, winEqL) # Glob poleward energy transport
LatMaxTransport = ('ntr', 'ntr', getLatMaxTransport, 'latmaxtransport', winNHemi, winNHemi) # latitude of the maximum energy transport

GlobSIC = ('sic', 'sic', getMean, 'globsic', winGlob, winGlob)
NHemiSIC = ('sic', 'sic', getMean, 'nhemisic', winNHemi, winNHemi) # Mean Sea Ice Cover

GlobDEP = ('gttd', 'var323', getMean, 'globdep', winGlob, winGlob) # Mean Diabatic Entropy Production
EqDEP = ('gttd', 'var323', getMean, 'eqdep', winEq, winEq) # Mean Diabatic Entropy Production at equator

MTG = ('tsa', 'tsa', getMerGradient, 'mtg', winNTrop, winNPole) # Zonal Mean Meridional Gradient of Surface Temperature

# Index choice
if indexChoice == 'globmst':
    indexDef = GlobMST
elif indexChoice == 'npolemst':
    indexDef = NPoleMST
elif indexChoice == 'ntropmst':
    indexDef = NTropMST
elif indexChoice == 'eqmst':
    indexDef = EqMST
elif indexChoice == 'npolentr':
    indexDef = NPoleNTR
elif indexChoice == 'eqlntr':
    indexDef = EqLNTR
elif indexChoice == 'energytransport':
    indexDef = EnergyTransport
elif indexChoice == 'latmaxtransport':
    indexDef = LatMaxTransport
elif indexChoice == 'globsic':
    indexDef = GlobSIC
elif indexChoice == 'nhemisic':
    indexDef = NHemiSIC
elif indexChoice == 'globdep':
    indexDef = GlobDEP
elif indexChoice == 'eqdep':
    indexDef = EqDEP
elif indexChoice == 'mtg':
    indexDef = MTG
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
    if firstYearOfFile[k] == 9901:
        print 'Processing last file beginning year %d' % firstYearOfFile[k]
        ncFilePath = '%s/%s_%s_%d_%05d_%05d.nc' \
            % (srcDir, indexDef[0], restartState, S,
               firstYearOfFile[k],
               firstYearOfFile[k] + yearsPerFile - 1 - 1)
    else:
        ncFilePath = '%s/%s_%s_%d_%05d_%05d.nc' \
            % (srcDir, indexDef[0], restartState, S,
               firstYearOfFile[k],
               firstYearOfFile[k] + yearsPerFile - 1)
    print 'Reading %s...' % ncFilePath
    dset = Dataset(ncFilePath, 'r')
    srcVar = dset.variables[indexDef[1]][:]

    # Process variable
    indexFile = indexDef[2](srcVar, lat, lon, WEIGHTS, indexDef[4], indexDef[5])
    if firstYearOfFile[k] == 9901:
        index[k*yearsPerFile*daysPerYear:((k+1)*yearsPerFile-1)*daysPerYear] = indexFile
    else:
        index[k*yearsPerFile*daysPerYear:(k+1)*yearsPerFile*daysPerYear] = indexFile

    # Close netCDF dataset
    dset.close()

# Save index
dstDir = '%s_%d/' % (restartState, S)
os.system('mkdir %s 2> /dev/null' % dstDir)
np.savetxt('%s/%s_%s_%d_%05d_%05d.txt' \
               % (dstDir, indexDef[3], restartState, S,
                  firstYear, lastYear), index)


