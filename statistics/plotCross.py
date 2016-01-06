import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import cm
import atmath

# Define the observable
srcDir = '../runPlasim/postprocessor/indices/'

SRng = np.array([1265, 1270, 1275, 1280, 1290, 1300, 1315, 1330, 1345, 1360])
restartStateRng = ['warm'] * SRng.shape[0]
lastYearRng = [9999] * SRng.shape[0]

firstYear = 101
yearsPerFile = 100
daysPerYear = 360
#indexChoice = ('nhemisic', 'energytransport',)
indexChoice = ('nhemisic', 'eqmst',)
ylabel = (r'NH SIC % and Eq MST',)
#indexChoice = ('eqmst', 'energytransport',)
#indexChoice = ('energytransport', 'latmaxtransport',)
#indexChoice = ('eqmst', 'latmaxtransport',)

# Case definition
spinupYears = 200  # Remove spinup period from time-series
spinup = spinupYears * daysPerYear
sampling = daysPerYear # (days^{-1})
sampFreq = 1 # (days^{-1})

# Plot settings
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
#lagMax = daysPerYear * 5
ccfRng = np.empty((SRng.shape[0], 2*lagMax+1))
PSDSlow = np.empty((SRng.shape[0], 600))
PSDSlowNormed = np.empty((SRng.shape[0], 600))

for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]

    # Read datasets
    resDir = '%s_%d/' % (restartState, S*10)
    dstDir = resDir
    indicesPath = '%s/%s/' % (srcDir, resDir)
    obsName = '%s_%d_%05d_%05d_anom' \
              % (restartState, S*10, firstYear, lastYearRng[k])
    indexFile1 = '%s_%s_%d_%05d_%05d.txt' % (indexChoice[0], restartState,
                                             S*10, firstYear, lastYearRng[k])
    indexFile2 = '%s_%s_%d_%05d_%05d.txt' % (indexChoice[1], restartState,
                                             S*10, firstYear, lastYearRng[k])
    os.system('mkdir stats %s %s/seasonal %s/anom_yearly 2> /dev/null' \
              % (dstDir, dstDir, dstDir))
    print 'Reading index file %s...' % indexFile1
    observable1 = np.loadtxt('%s/%s' % (indicesPath, indexFile1))
    print 'Reading index file %s...' % indexFile2
    observable2 = np.loadtxt('%s/%s' % (indicesPath, indexFile2))
    obsName += '_%s' % indexChoice[0] +'_%s' % indexChoice[1]
    ntFull = observable1.shape[0]

    # Get time steps array
    observable1 = observable1[spinup:]
    observable1Yearly = np.convolve(observable1, np.ones((sampling,)) \
                                    / sampling)[(sampling-1)::sampling]
    nt = ntFull - spinup
    observable2 = observable2[spinup:]
    observable2Yearly = np.convolve(observable2, np.ones((sampling,)) \
                                    / sampling)[(sampling-1)::sampling]

    # Get cross-correlation function
    ccfRng[k] = atmath.ccf(observable1Yearly, observable2Yearly, lagMax=lagMax)
    np.savetxt('%s/stats/ccf_yearly_%s.txt' % (dstDir, obsName), ccfRng[k])
                
#     # Get cross-periodogram
#     print 'Getting cross-periodogram...'
#     window = np.hamming(nt)
#     # Get nearest larger power of 2
#     if np.log2(nt) != int(np.log2(nt)):
#         nfft = 2**(int(np.log2(nt)) + 1)
#     else:
#         nfft = nt
#     # Get frequencies and shift zero frequency to center
#     freq = np.fft.fftfreq(nfft, d=1./sampFreq)
#     freq = np.fft.fftshift(freq)
#     freqYear = freq * daysPerYear

#     # Apply window and remove sample mean
#     ts1Windowed = observable1 * window
#     ts1Windowed -= ts1Windowed.mean()
#     ts2Windowed = observable2 * window
#     ts2Windowed -= ts2Windowed.mean()
#     # Fourier transform and shift zero frequency to center
#     fts1 = np.fft.fft(ts1Windowed, nfft, 0)
#     fts1 = np.fft.fftshift(fts1)
#     fts2 = np.fft.fft(ts2Windowed, nfft, 0)
#     fts2 = np.fft.fftshift(fts2)
#     # Get dimensional periodogram (summing to the variance)
#     perio = np.abs(fts1 * fts2) / nt / nfft
#     # Get spectral density (normalized by the variance)
#     perioSum = perio.sum()
#     perio /= perioSum
        
#     # Plot
#  #   perio[perio < 10**(-5)] = 10**(-5)
#     print 'Plotting...'
#     swin = '_log'
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)
#     ax.plot(freqYear[nfft/2+1:], perio[nfft/2+1:], '.', markersize=5)
#     ax.set_xscale('log')
#     ax.set_yscale('log')
# #    ax.set_xlim(freqYear[freq.shape[0]/2+4000], freqYear[freq.shape[0]/2+420000])
#     ax.set_ylim(10**(-12), 10**(-1))
#     ax.set_xlabel(r'years$^{-1}$', fontsize=fs_latex)
#     ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
#     plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
#     plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
#     plt.title('Periodogram  for case %s_%d' % (restartState, S))
#     fig.savefig('%s/perio/crossperio%s_%s.%s' % (dstDir, swin, obsName, figFormat),
#                 bbox_inches='tight', dpi=dpi)


# fig = plt.figure()
# ax = fig.add_subplot(111)
# for k in np.arange(SRng.shape[0]/2):
#     S = SRng[k]
#     ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlow[k], label=str(S), linestyle='-')
# for k in np.arange(SRng.shape[0]/2, SRng.shape[0]):
#     S = SRng[k]
#     ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlow[k], label=str(S), linestyle='--')
# ax.legend(loc='lower right')
# fig.savefig('stats/crossPSDSlow_%s.%s' % (indexChoice[0], figFormat),
#             bbox_inches='tight', dpi=dpi)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# for k in np.arange(SRng.shape[0]/2):
#     S = SRng[k]
#     ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlowNormed[k],
#             label=str(S), linestyle='-')
# for k in np.arange(SRng.shape[0]/2, SRng.shape[0]):
#     S = SRng[k]
#     ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlowNormed[k],
#             label=str(S), linestyle='--')
# ax.legend(loc='lower right')
# fig.savefig('stats/crossPSDSlowNormed_%s.%s' % (indexChoice[0], figFormat),
#             bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
for k in np.arange(0, SRng.shape[0], 2):
    S = SRng[k]
    ax.plot(np.arange(-lagMax, lagMax+1), ccfRng[k], label=str(S),
            linestyle='-', linewidth=2)
    S = SRng[k+1]
    ax.plot(np.arange(-lagMax, lagMax+1), ccfRng[k+1], label=str(S),
            linestyle='--', linewidth=2)
ax.set_xlim(-lagMax, lagMax)
ax.set_xlabel('Lag (years)', fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
#ax.set_ylim(-0.05, 1.)
ax.set_ylabel(r'Cross-correlation of %s' % ylabel[0], fontsize=fs_latex)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
ax.legend(loc='upper right', title=r'$S$ $(W \cdot m^{-2})$')
fig.savefig('interComparison/ccf_%s_%s.%s' % (indexChoice[0], indexChoice[1], figFormat),
            bbox_inches='tight', dpi=dpi)

