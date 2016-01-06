import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import cm
import atmath

# Define the observable
srcDir = '../runPlasim/postprocessor/indices/'

# SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433,
#                  1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
# restartStateRng = np.concatenate((['cold']*8, ['warm']*8), 0)

SRng = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartStateRng = ['warm']*8

#SRng = np.array([1263, 1265])
#restartStateRng = ['warm']*2


firstYear = 101
lastYear = 4200
yearsPerFile = 100
daysPerYear = 360
#indexChoice = ('globmst',)
#indexChoice = ('npolemst',)
#indexChoice = ('globdep',)
#indexChoice = ('eqdep',)
#indexChoice = ('MTG',)
#indexChoice = ('areabelowtf20nhemi',)
indexChoice = ('areabelowtfnhemi',)

# Case definition
spinupYears = 100  # Remove spinup period from time-series
spinup = spinupYears * daysPerYear
sampFreq = 1 # (days^{-1})

# Plot settings
fs_default = 'x-large'
fs_latex = 'xx-large'
fs_xlabel = fs_default
fs_ylabel = fs_default
fs_xticklabels = fs_default
fs_yticklabels = fs_default
fs_legend_title = fs_default
fs_legend_labels = fs_default
fs_cbar_label = fs_default
#            figFormat = 'eps'
figFormat = 'png'
dpi = 300

varRng = np.empty((SRng.shape[0],))
skewRng = np.empty((SRng.shape[0],))
kurtRng = np.empty((SRng.shape[0],))
lagMax = 80
#lagMax = daysPerYear * 5
ccfRng = np.empty((SRng.shape[0], lagMax*2+1))

for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]

    # Create directories
    resDir = '%s_%s/' % (restartState, S)
    dstDir = resDir
    indicesPath = '%s/%s/' % (srcDir, resDir)
    os.system('mkdir stats %s %s/seasonal %s/anom 2> /dev/null' % (dstDir, dstDir, dstDir))
    
    # Read datasets
    obsName = '%s_%d_%05d_%05d_anom' % (restartState, S, firstYear, lastYear)
    indexFile = '%s_%s_%d_%05d_%05d.txt' \
                % (indexChoice[0], restartState, S, firstYear, lastYear)
    print 'Reading index file %s...' % indexFile
    observable = np.loadtxt('%s/%s' % (indicesPath, indexFile))
    ntFull = observable.shape[0]
    obsName += '_%s' % indexChoice[0]

    # Get time steps array
    time = np.arange(spinup, ntFull)
    nt = ntFull - spinup
    observable = observable[spinup:]
    seasonal = np.empty((daysPerYear,))
    anom = np.empty((nt,))
    for day in np.arange(daysPerYear):
        seasonal[day] = observable[day::daysPerYear].mean()
        anom[day::daysPerYear] = observable[day::daysPerYear] - seasonal[day]

    varRng[k] = anom.var()
    skewRng[k] = stats.skew(anom)
    kurtRng[k] = stats.kurtosis(anom)
    ccfRng[k] = atmath.ccf(anom, anom, lagMax=lagMax)
        
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(1, daysPerYear+1), seasonal)
    ax.set_xlabel(r'days', fontsize=fs_latex)
    ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    plt.title('Seasonal cycle for case %s_%d\n\sigma = %.5f' % (restartState, S, seasonal.std()))
    fig.savefig('%s/seasonal/seasonal_%s.%s' % (dstDir, obsName, figFormat),
                bbox_inches='tight', dpi=dpi)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(time[200*daysPerYear:203*daysPerYear], anom[200*daysPerYear:203*daysPerYear])
    ax.set_xlabel(r'days', fontsize=fs_latex)
    ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    plt.title('Anomalies for case %s_%d\n\sigma = %.5f' % (restartState, S, anom.std()))
    fig.savefig('%s/anom/anom_%s.%s' % (dstDir, obsName, figFormat),
                bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, varRng)
fig.savefig('stats/variance_%s.%s' % (indexChoice[0], figFormat), bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, skewRng)
fig.savefig('stats/skewness_%s.%s' % (indexChoice[0], figFormat), bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, kurtRng)
fig.savefig('stats/kurtosis_%s.%s' % (indexChoice[0], figFormat), bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
for k in np.arange(SRng.shape[0]/2):
    S = SRng[k]
    ax.plot(np.arange(-lagMax, lagMax+1), ccfRng[k], label=str(S), linestyle='-')
for k in np.arange(SRng.shape[0]/2, SRng.shape[0]):
    S = SRng[k]
    ax.plot(np.arange(-lagMax, lagMax+1), ccfRng[k], label=str(S), linestyle='--')
ax.legend(loc='upper right')
ax.set_xlim(0, lagMax)
ax.set_ylim(-0.05, 1.)
fig.savefig('stats/acf_%s.%s' % (indexChoice[0], figFormat), bbox_inches='tight', dpi=dpi)
