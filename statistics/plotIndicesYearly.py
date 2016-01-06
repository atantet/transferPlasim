import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import atmath

# Define the observable
srcDir = '../runPlasim/postprocessor/indices/'

# SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433,
#                  1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
# restartStateRng = np.concatenate((['cold']*8, ['warm']*8), 0)

SRng = np.array([1265, 1270, 1275, 1280, 1290, 1300, 1315, 1330, 1345, 1360])
restartStateRng = ['warm'] * SRng.shape[0]
lastYearRng = [9999] * SRng.shape[0]

# SRng1 = np.array([1263, 1264, 1265, 1270, 1275, 1280, 1290, 1300, 1315, 1330, 1345, 1360, 1400, 1435])
# restartStateRng1 = ['warm'] * SRng1.shape[0]
# lastYearRng1 = [9999] * SRng1.shape[0]
# SRng2 = np.array([1262.8, 1262.9, 1263.5, 1264.5])
# restartStateRng2 = ['warm'] * SRng2.shape[0]
# lastYearRng2 = [1500] * SRng2.shape[0]
# SRng = np.concatenate((SRng1, SRng2), 0)
# restartStateRng = np.concatenate((restartStateRng1, restartStateRng2), 0)
# lastYearRng = np.concatenate((lastYearRng1, lastYearRng2), 0)
# isort = np.argsort(SRng)
# SRng = SRng[isort]
# restartStateRng = restartStateRng[isort]
# lastYearRng = lastYearRng[isort]

# SRng1 = np.array([1263, 1264, 1265, 1270, 1275, 1280])
# restartStateRng1 = ['warm'] * SRng1.shape[0]
# lastYearRng1 = [9999] * SRng1.shape[0]
# SRng2 = np.array([1262.8, 1262.9, 1263.5, 1264.5])
# restartStateRng2 = ['warm'] * SRng2.shape[0]
# lastYearRng2 = [6300] * SRng2.shape[0]
# SRng = np.concatenate((SRng1, SRng2), 0)
# restartStateRng = np.concatenate((restartStateRng1, restartStateRng2), 0)
# lastYearRng = np.concatenate((lastYearRng1, lastYearRng2), 0)
# isort = np.argsort(SRng)
# SRng = SRng[isort]
# restartStateRng = restartStateRng[isort]
# lastYearRng = lastYearRng[isort]

# SRng1 = np.array([1265, 1270, 1275, 1280, 1290, 1300, 1315, 1330, 1345, 1360])
# restartStateRng1 = ['warm'] * SRng1.shape[0]
# lastYearRng1 = [9999] * SRng1.shape[0]
# SRng2 = np.array([1267.5, 1272.5, 1285, 1295])
# restartStateRng2 = ['warm'] * SRng2.shape[0]
# lastYearRng2 = [4000] * SRng2.shape[0]
# SRng = np.concatenate((SRng1, SRng2), 0)
# restartStateRng = np.concatenate((restartStateRng1, restartStateRng2), 0)
# lastYearRng = np.concatenate((lastYearRng1, lastYearRng2), 0)
# isort = np.argsort(SRng)
# SRng = SRng[isort]
# restartStateRng = restartStateRng[isort]
# lastYearRng = lastYearRng[isort]

firstYear = 101
yearsPerFile = 100
daysPerYear = 360
#indexChoice = ('globmst',)
#indexChoice = ('npolemst',)
#indexChoice = ('eqmst',)
#ylabel = (r'Eq MST', 'K')
#indexChoice = ('MTG',)
#ylabel = (r'MTD', 'K')
indexChoice = ('energytransport',)
ylabel = (r'GET', 'PW')
#indexChoice = ('energybalance',)
#ylabel = (r'GEB', 'PW')
#indexChoice = ('latmaxtransport',)
#ylabel = (r'LMET', 'degree')
#indexChoice = ('eqlntr',)
#indexChoice = ('globsic',)
#indexChoice = ('nhemisic',)
#ylabel = (r'NH SIC %', '')
#indexChoice = ('globdep',)
#indexChoice = ('eqdep',)

# Case definition
spinupYears = 100  # Remove spinup period from time-series
spinup = spinupYears * daysPerYear
sampling = daysPerYear # (days^{-1})

# Plot settings
xlabel = r'$S$ $(W \cdot m^{-2})$' 
xticklabels = SRng.astype(str)
xticklabels[[0, 2]] = ''
xlim = (1262.5, 1362.5)
ls = '--'
lw = 2
marker = 'o'
msize = 10
fs_default = 'x-large'
fs_latex = 'xx-large'
fs_xlabel = fs_default
fs_ylabel = fs_default
fs_xticklabels = 'large'
fs_yticklabels = 'large'
fs_legend_title = fs_default
fs_legend_labels = fs_default
fs_cbar_label = fs_default
#            figFormat = 'eps'
figFormat = 'png'
dpi = 300

meanRng = np.empty((SRng.shape[0],))
varRng = np.empty((SRng.shape[0],))
skewRng = np.empty((SRng.shape[0],))
kurtRng = np.empty((SRng.shape[0],))
lagMax = 100
#lagMax = daysPerYear * 5
ccfRng = np.empty((SRng.shape[0], lagMax+1))

os.system('mkdir interComparison 2> /dev/null')
for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]

    # Read datasets
    resDir = '%s_%d/' % (restartState, S*10)
    dstDir = resDir
    indicesPath = '%s/%s/' % (srcDir, resDir)
    obsName = '%s_%d_%05d_%05d_anom' \
              % (restartState, S*10, firstYear, lastYearRng[k])
    indexFile = '%s_%s_%d_%05d_%05d.txt' % (indexChoice[0], restartState,
                                            S*10, firstYear, lastYearRng[k])
    os.system('mkdir %s %s/anom_yearly %s/stats/ 2> /dev/null' \
              % (dstDir, dstDir, dstDir))
    print 'Reading index file %s...' % indexFile
    observable = np.loadtxt('%s/%s' % (indicesPath, indexFile))
    
    ntFull = observable.shape[0]
    obsName += '_%s' % indexChoice[0]

    # Get time steps array
    observable = observable[spinup:]
    observable = np.convolve(observable,
                             np.ones((sampling,))/sampling)[(sampling-1)::sampling]
    nt = observable.shape[0]
    time = np.arange(nt)

    meanRng[k] = observable.mean()
    varRng[k] = observable.var()
    skewRng[k] = stats.skew(observable)
    kurtRng[k] = stats.kurtosis(observable)
    ccfRng[k] = atmath.ccf(observable, observable, lagMax=lagMax)[lagMax:]
    np.savetxt('%s/stats/mean_yearly_%s.txt' % (dstDir, obsName),
               (meanRng[k],))
    np.savetxt('%s/stats/var_yearly_%s.txt' % (dstDir, obsName),
               (varRng[k],))
    np.savetxt('%s/stats/skew_yearly_%s.txt' % (dstDir, obsName),
               (skewRng[k],))
    np.savetxt('%s/stats/kurt_yearly_%s.txt' % (dstDir, obsName),
               (kurtRng[k],))
    np.savetxt('%s/stats/acf_yearly_%s.txt' % (dstDir, obsName),
               ccfRng[k])
                
    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1)
    # ax.plot(time[30:230], observable[30:230])
    # ax.set_xlabel('years', fontsize=fs_latex)
    # ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
    # plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    # plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    # plt.title('Observable for case %s_%d\n$\sigma = %.5f$' \
    #           % (restartState, S, observable.std()))
    # fig.savefig('%s/anom_yearly/yearly_%s.%s' % (dstDir, obsName, figFormat),
    #             bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, meanRng, linestyle=ls, marker=marker, markersize=msize)
ax.set_xlabel(xlabel, fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
ax.set_xlim(xlim)
ax.set_xticks(SRng)
ax.set_xticklabels(xticklabels)
if ylabel[1] !=  '':
    ax.set_ylabel(r'Mean of %s ($%s$)' % ylabel, fontsize=fs_latex)
else:
    ax.set_ylabel(r'Mean of %s' % ylabel[0], fontsize=fs_latex)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
fig.savefig('interComparison/mean_yearly_%s.%s' % (indexChoice[0], figFormat),
            bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, varRng, linestyle=ls, marker=marker, markersize=msize)
ax.set_xlabel(xlabel, fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
ax.set_xlim(xlim)
ax.set_xticks(SRng)
ax.set_xticklabels(xticklabels)
if ylabel[1] !=  '':
    ax.set_ylabel(r'Variance of %s ($%s^2$)' % ylabel, fontsize=fs_latex)
else:
    ax.set_ylabel(r'Variance of %s' % ylabel[0], fontsize=fs_latex)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
fig.savefig('interComparison/variance_yearly_%s.%s' % (indexChoice[0], figFormat),
            bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, skewRng, linestyle=ls, marker=marker, markersize=msize)
ax.set_xlabel(xlabel, fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
ax.set_xlim(xlim)
ax.set_xticks(SRng)
ax.set_xticklabels(xticklabels)
ax.set_ylabel(r'Skewness of %s' % ylabel[0], fontsize=fs_latex)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
fig.savefig('interComparison/skewness_yearly_%s.%s' % (indexChoice[0], figFormat),
            bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, kurtRng, linestyle=ls, marker=marker, markersize=msize)
ax.set_xlabel(xlabel, fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
ax.set_xlim(xlim)
ax.set_xticks(SRng)
ax.set_xticklabels(xticklabels)
ax.set_ylabel(r'Kurtosis of %s' % ylabel[0], fontsize=fs_latex)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
fig.savefig('interComparison/kurtosis_yearly_%s.%s' % (indexChoice[0], figFormat),
            bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
for k in np.arange(0, SRng.shape[0], 2):
    S = SRng[k]
    ax.plot(np.arange(0., lagMax+1), ccfRng[k], label=str(S),
            linestyle='-', linewidth=2)
    S = SRng[k+1]
    ax.plot(np.arange(0., lagMax+1), ccfRng[k+1], label=str(S),
            linestyle='--', linewidth=2)
ax.set_xlim(0, lagMax)
ax.set_xlabel('Lag (years)', fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
ax.set_ylim(-0.05, 1.)
ax.set_ylabel(r'Auto-correlation of %s' % ylabel[0], fontsize=fs_latex)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
ax.legend(loc='upper right', title=r'$S$ $(W \cdot m^{-2})$')
fig.savefig('interComparison/acf_%s.%s' % (indexChoice[0], figFormat),
            bbox_inches='tight', dpi=dpi)

