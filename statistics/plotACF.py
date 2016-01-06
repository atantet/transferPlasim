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

firstYear = 101
yearsPerFile = 100
daysPerYear = 360
#indexChoice = ('globmst',)
#indexChoice = ('npolemst',)
#indexChoice = ('eqmst',)
#ylabel = (r'Eq MST', 'K')
#indexChoice = ('MTG',)
#ylabel = (r'MTD', 'K')
#indexChoice = ('energytransport',)
#ylabel = (r'GET', 'PW')
#indexChoice = ('latmaxtransport',)
#ylabel = (r'LMET', 'degree')
#indexChoice = ('eqlntr',)
#indexChoice = ('globsic',)
indexChoice = ('nhemisic',)
ylabel = (r'NH SIC %', '')
#indexChoice = ('globdep',)
#indexChoice = ('eqdep',)

# Case definition
spinupYears = 100  # Remove spinup period from time-series
spinup = spinupYears * daysPerYear
sampling = daysPerYear # (days^{-1})

# Plot settings
lagMaxPlot = 50
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

lagMax = 100
ccfRng = np.empty((SRng.shape[0], lagMax+1))

os.system('mkdir interComparison 2> /dev/null')
for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]

    # Read datasets
    resDir = '%s_%d/' % (restartState, S*10)
    dstDir = resDir
    indicesPath = '%s/%s/' % (srcDir, resDir)
    obsName = '%s_%d_%05d_%05d_anom' % (restartState, S*10, firstYear, lastYearRng[k])
    indexFile = '%s_%s_%d_%05d_%05d.txt' \
                % (indexChoice[0], restartState, S*10, firstYear, lastYearRng[k])
    obsName += '_%s' % indexChoice[0]
    ccfRng[k] = np.loadtxt('%s/stats/acf_yearly_%s.txt' % (dstDir, obsName))
                
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

fig = plt.figure()
ax = fig.add_subplot(111)
for k in np.arange(0, SRng.shape[0], 2):
    S = SRng[k]
    ax.plot(np.arange(0., lagMaxPlot+1), ccfRng[k][:lagMaxPlot+1], label=str(S),
            linestyle='-', linewidth=2)
    S = SRng[k+1]
    ax.plot(np.arange(0., lagMaxPlot+1), ccfRng[k+1][:lagMaxPlot+1], label=str(S),
            linestyle='--', linewidth=2)
ax.set_xlim(0, lagMaxPlot)
ax.set_xlabel('Lag (years)', fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
ax.set_ylim(1.e-2, 1.)
ax.set_yscale('log')
ax.set_ylabel(r'Auto-correlation of %s' % ylabel[0], fontsize=fs_latex)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
ax.legend(loc='upper right', title=r'$S$ $(W \cdot m^{-2})$')
fig.savefig('interComparison/acfLog_%s.%s' % (indexChoice[0], figFormat),
            bbox_inches='tight', dpi=dpi)

