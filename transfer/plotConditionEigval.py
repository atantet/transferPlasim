import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
import atmath

    # Case definition
SRng = np.array([1265, 1270, 1275, 1280, 1300, 1315, 1330, 1345, 1360])
restartState = "warm"
lastYear = 9999
firstYear = 101
processing = "_yearly"
indexChoice = ["nhemisic"]
dim = len(indexChoice)

# Grid definition
nx = 50
nSTD = 5

# Lags
nLags = 1
tauDimRng = np.array([1])

nev = 10

# Plot
fs_default = 'x-large'
fs_latex = 'xx-large'
fs_xlabel = fs_default
fs_ylabel = fs_default
fs_xticklabels = fs_default
fs_yticklabels = fs_default
fs_legend_title = fs_default
fs_legend_labels = fs_default
fs_cbar_label = fs_default
msize = 48
scattersize = 12
#scattersize = 36
#            figFormat = 'eps'
figFormat = 'png'
dpi = 300
readMap = False


conditionSecondEigval = np.empty((SRng.shape[0],))
for iS in np.arange(SRng.shape[0]):
    S = SRng[iS]
    resDir = "%s_%d" % (restartState, np.round(S*10))
    dstDir = resDir
    caseName = "%s_%05d_%05d" % (resDir, firstYear, lastYear)
    obsName = caseName
    gridPostfix = "N"
    N = 1
    for d in np.arange(dim):
        N *= nx
        obsName += ("_%s" % indexChoice[d])
        if d > 0:
            gridPostfix += ("x%d" % nx)
        else:
            gridPostfix += ("%d" % nx)
    gridPostfix = "%s_%s_%s_%dstd" % (processing, obsName, gridPostfix, nSTD)

    for lag in np.arange(tauDimRng.shape[0]):
        tauDim = tauDimRng[lag]
        postfix = "%s_tau%02d" % (gridPostfix, tauDim)

    # Get condition number
    print 'Getting condition number...'
    conditionSecondEigval[iS] = np.loadtxt('%s/spectrum/condition/conditionEigval_nev%d%s.txt' % (dstDir, nev, postfix))[1]
        
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(SRng, conditionSecondEigval)
fig.savefig('conditionSecondeigval_nev%d%s.%s' \
            % (nev, postfix, figFormat),
            bbox_inches='tight', dpi=dpi)
