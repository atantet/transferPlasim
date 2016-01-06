import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define the observable
srcDir = '../runPlasim/postprocessor/indices/'

# SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433,
#                  1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
# restartStateRng = np.concatenate((['cold']*8, ['warm']*8), 0)

SRng = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartStateRng = ['warm']*8

SRng1 = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartStateRng1 = ['warm'] * SRng1.shape[0]
lastYearRng1 = [9999] * SRng1.shape[0]
SRng2 = np.array([1264, 1275, 1290, 1315, 1345, 1400])
restartStateRng2 = ['warm'] * SRng2.shape[0]
lastYearRng2 = [5500] * SRng2.shape[0]
SRng = np.concatenate((SRng1, SRng2), 0)
restartStateRng = np.concatenate((restartStateRng1, restartStateRng2), 0)
lastYearRng = np.concatenate((lastYearRng1, lastYearRng2), 0)
isort = np.argsort(SRng)
SRng = SRng[isort]
restartStateRng = restartStateRng[isort]
lastYearRng = lastYearRng[isort]

firstYear = 101
yearsPerFile = 100
daysPerYear = 360
#indexChoice = ('globmst',)
#indexChoice = ('npolemst',)
#indexChoice = ('eqmst',)
#indexChoice = ('eqlntr',)
indexChoice = ('nhemisic',)
#indexChoice = ('globdep',)
#indexChoice = ('eqdep',)
#indexChoice = ('MTG',)
#indexChoice = ('areabelowtf20nhemi',)

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

PSDSlow = np.empty((SRng.shape[0], 600))
PSDSlowNormed = np.empty((SRng.shape[0], 600))

for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]

    # Read datasets
    obsName = '%s_%d_%05d_%05d_anom' % (restartState, S, firstYear, lastYearRng[k])
    indexFile = '%s_%s_%d_%05d_%05d.txt' \
                % (indexChoice[0], restartState, S, firstYear, lastYearRng[k])
    resDir = '%s_%d/' % (restartState, S)
    dstDir = resDir
    indicesPath = '%s/%s/' % (srcDir, resDir)
    if os.path.isfile('%s/%s' % (indicesPath, indexFile)):
        os.system('mkdir %s %s/perio 2> /dev/null' % (dstDir, dstDir))
        print 'Reading index file %s...' % indexFile
        observable = np.loadtxt('%s/%s' % (indicesPath, indexFile))
    else:
        resDir = '%s_%d/' % (restartState, S*10)
        dstDir = resDir
        indicesPath = '%s/%s/' % (srcDir, resDir)
        obsName = '%s_%d_%05d_%05d_anom' % (restartState, S*10, firstYear, lastYearRng[k])
        indexFile = '%s_%s_%d_%05d_%05d.txt' \
                    % (indexChoice[0], restartState, S*10, firstYear, lastYearRng[k])
        os.system('mkdir %s %s/perio 2> /dev/null' % (dstDir, dstDir))
        print 'Reading index file %s...' % indexFile
        observable = np.loadtxt('%s/%s' % (indicesPath, indexFile))
    ntFull = observable.shape[0]
    obsName += '_%s' % indexChoice[0]

    # Get time steps array
    time = np.arange(spinup, ntFull)
    nt = ntFull - spinup
    observable = observable[spinup:]

    # Get periodogram
    print 'Getting periodogram...'
    window = np.hamming(nt)
    # Get nearest larger power of 2
    if np.log2(nt) != int(np.log2(nt)):
        nfft = 2**(int(np.log2(nt)) + 1)
    else:
        nfft = nt
    # Get frequencies and shift zero frequency to center
    freq = np.fft.fftfreq(nfft, d=1./sampFreq)
    freq = np.fft.fftshift(freq)
    freqYear = freq * daysPerYear

    # Apply window and remove sample mean
    tsWindowed = observable * window
    tsWindowed -= tsWindowed.mean()
    # Fourier transform and shift zero frequency to center
    fts = np.fft.fft(tsWindowed, nfft, 0)
    fts = np.fft.fftshift(fts)
    # Get dimensional periodogram (summing to the variance)
    perio = np.abs(fts)**2 / nt / nfft
    # Get spectral density (normalized by the variance)
    perioSum = perio.sum()
    perio /= perioSum
        
    # Cumulative periodogram, normalized
    PSDSlow[k] = perio[nfft/2:nfft/2+PSDSlow.shape[1]].cumsum() * perioSum
    PSDSlowNormed[k] = perio[nfft/2:nfft/2+PSDSlow.shape[1]].cumsum()
#    print 'Energy in periods smaller than two years:', PSDSlow[k]
    
#   # Plot
 #   perio[perio < 10**(-5)] = 10**(-5)
    print 'Plotting...'
    swin = '_log'
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(freqYear[nfft/2+1:], perio[nfft/2+1:], '.', markersize=5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-4, 2e2)
    ax.set_ylim(10**(-12), 10**(-1))
    ax.set_xlabel(r'years$^{-1}$', fontsize=fs_latex)
    ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    plt.title('Periodogram  for case %s_%d' % (restartState, S))
    fig.savefig('%s/perio/perio%s_%s.%s' % (dstDir, swin, obsName, figFormat),
                bbox_inches='tight', dpi=dpi)

#    swin = '_decacdal'
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)
#     ax.plot(freqYear[nfft/2+20:nfft/2+6000], np.log10(perio)[nfft/2+20:nfft/2+6000],
#             '.', markersize=5)
#     #    ax.set_xscale('log')
#     #    ax.set_yscale('log')
#  #   ax.set_xlim(freqYear[nfft/2+20], freqYear[nfft/2+4000])
#     ax.set_ylim(-12, 0)
#     ax.set_xlabel(r'years$^{-1}$', fontsize=fs_latex)
#     ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
#     plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
#     plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
#     plt.title('Periodogram  for case %s_%d' % (restartState, S))
#     fig.savefig('%s/perio/perio%s_%s.%s' % (dstDir, swin, obsName, figFormat),
#                 bbox_inches='tight', dpi=dpi)

#     swin = '_seasonal'
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)
#     ax.plot(freq[nfft/2+4000:nfft/2+100000], np.log10(perio)[nfft/2+4000:nfft/2+100000],
#             '.', markersize=5)
#     #    ax.set_xscale('log')
#     #    ax.set_yscale('log')
# #    ax.set_xlim(freq[freq.shape[0]/2+4000], freq[freq.shape[0]/2+420000])
#     ax.set_ylim(-12, 0)
#     ax.set_xlabel(r'days$^{-1}$', fontsize=fs_latex)
#     ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
#     plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
#     plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
#     plt.title('Periodogram  for case %s_%d' % (restartState, S))
#     fig.savefig('%s/perio/perio%s_%s.%s' % (dstDir, swin, obsName, figFormat),
#                 bbox_inches='tight', dpi=dpi)


#     swin = '_synoptic'
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)
#     ax.plot(freq[freq.shape[0]/2+100000:],
#             np.log10(perio)[freq.shape[0]/2+100000:], '.', markersize=5)
#     #    ax.set_xscale('log')
#     #    ax.set_yscale('log')
#     #    ax.set_xlim(freq[freq.shape[0]/2+4000], freq[freq.shape[0]/2+420000])
#     ax.set_ylim(-12, 0)
#     ax.set_xlabel(r'days$^{-1}$', fontsize=fs_latex)
#     ax.set_ylabel(indexChoice[0], fontsize=fs_latex)
#     plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
#     plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
#     plt.title('Periodogram  for case %s_%d' % (restartState, S))
#     fig.savefig('%s/perio/perio%s_%s.%s' % (dstDir, swin, obsName, figFormat),
#                 bbox_inches='tight', dpi=dpi)

# Plot energy in low frequencies
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(SRng, PSDSlow, linewidth=2)
# fig.savefig('stats/PSDSlow_%s.%s' % (indexChoice[0], figFormat),
#             bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
for k in np.arange(SRng.shape[0]/2):
    S = SRng[k]
    ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlow[k], label=str(S), linestyle='-')
for k in np.arange(SRng.shape[0]/2, SRng.shape[0]):
    S = SRng[k]
    ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlow[k], label=str(S), linestyle='--')
ax.legend(loc='lower right')
fig.savefig('stats/PSDSlow_%s.%s' % (indexChoice[0], figFormat), bbox_inches='tight', dpi=dpi)

fig = plt.figure()
ax = fig.add_subplot(111)
for k in np.arange(SRng.shape[0]/2):
    S = SRng[k]
    ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlowNormed[k],
            label=str(S), linestyle='-')
for k in np.arange(SRng.shape[0]/2, SRng.shape[0]):
    S = SRng[k]
    ax.plot(freqYear[nfft/2:nfft/2+PSDSlow.shape[1]], PSDSlowNormed[k],
            label=str(S), linestyle='--')
ax.legend(loc='lower right')
fig.savefig('stats/PSDSlowNormed_%s.%s' % (indexChoice[0], figFormat),
            bbox_inches='tight', dpi=dpi)

