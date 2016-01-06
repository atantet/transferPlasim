import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
import atmath

# Case definition
S = 1265
restartState = "warm"
lastYear = 9999
firstYear = 101
processing = "_yearly"
indexChoice = ["nhemisic", "eqmst"]
dim = len(indexChoice)

# Grid definition
nx = 25
nSTD = 5

# Lags
nLags = 1
tauDimRng = np.array([1])

nev = 25
nevPlot = 3
#plotAdjoint = False
plotAdjoint = True

resDir = "%s_%d" % (restartState, np.round(S*10))
dstDir = resDir
caseName = "%s_%05d_%05d" % (resDir, firstYear, lastYear)
os.system('mkdir %s/ccf %s/spectrum/eigval/figs %s/spectrum/eigvec/figs 2> /dev/null' \
          % (dstDir, dstDir, dstDir))

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


# Read grid
gridFile = '%s/grid/grid%s.txt' % (dstDir, gridPostfix)
f = open(gridFile, 'r')
bounds = []
coord = []
for k in np.arange(dim):
    bounds.append(np.array(f.readline().split()).astype(float))
    coord.append((bounds[k][1:] + bounds[k][:-1]) / 2)
f.close()
X, Y = np.meshgrid(coord[0], coord[1])

# Plot
levels = 20
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
gridXlim = [coord[0].min(), coord[0].max()]
gridYlim = [coord[1].min(), coord[1].max()]


for lag in np.arange(tauDimRng.shape[0]):
    tauDim = tauDimRng[lag]
    postfix = "%s_tau%02d" % (gridPostfix, tauDim)

    print 'Readig spectrum...'
    EigValFile = '%s/spectrum/eigval/eigval_nev%d%s.txt' % (dstDir, nev, postfix)
    EigVecFile = '%s/spectrum/eigvec/eigvec_nev%d%s.txt' % (dstDir, nev, postfix)
    eigval = np.loadtxt(EigValFile)
    eigval = eigval[:, 0] + eigval[:, 1]*1j
    eigvec = np.loadtxt(EigVecFile)
    eigvec = eigvec[::2] + eigvec[1::2]*1j
    isort = np.argsort(np.abs(eigval))[::-1]
    eigval = eigval[isort]
    eigvec = eigvec[isort]

    if plotAdjoint:
        EigValAdjointFile = '%s/spectrum/eigval/eigvalAdjoint_nev%d%s.txt' % (dstDir, nev, postfix)
        EigVecAdjointFile = '%s/spectrum/eigvec/eigvecAdjoint_nev%d%s.txt' % (dstDir, nev, postfix)
        eigvalAdjoint = np.loadtxt(EigValAdjointFile)
        eigvalAdjoint = eigvalAdjoint[:, 0] + eigvalAdjoint[:, 1]*1j
        eigvecAdjoint = np.loadtxt(EigVecAdjointFile)
        # From the transpose we get the conjugate of the adjoint eigenvectors
        # so we take back the conjugate
        eigvecAdjoint = eigvecAdjoint[::2] - eigvecAdjoint[1::2]*1j
        isort = np.argsort(np.abs(eigvalAdjoint))[::-1]
        eigvalAdjoint = eigvalAdjoint[isort]
        eigvecAdjoint = eigvecAdjoint[isort]
    nevSingle = eigval.shape[0]
        
    eigvalGen = np.empty((nevSingle,), dtype=complex)
    ev = 0
    for count in np.arange(eigval.shape[0]):
        eigvalGen[ev] = (np.log(np.abs(eigval[count])) \
                         + np.angle(eigval[count]) * 1j) / tauDim
        ev += 1
        if ev >= nevSingle:
            break
        if eigval[count].imag != 0:
            eigvalGen[ev] = np.conjugate(eigvalGen[ev-1])
            ev +=1
            if ev >= nevSingle:
                break


    # Plot spectrum
    print 'Plotting spectrum slowest rate ', -1. / eigvalGen[1].real
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(eigvalGen[1:].real, eigvalGen[1:].imag, c='k', s=msize, marker='o')
    ax.scatter(eigvalGen[0].real, eigvalGen[0].imag, c='r', s=msize, marker='o')
    ax.set_xlabel(r'$\Re(\zeta_i)$', fontsize=fs_latex)
    ax.set_ylabel(r'$\Im(\zeta_i)$', fontsize=fs_latex)
    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    #ax.set_title('%d-time-step spectrum for %s\nSlowest time-scale: %.1f' \
        #    % (tau, srcPostfix, -1. / rate[0]))
    ax.set_xlim(-1.2, 0.02)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.text(xlim[1] - (xlim[1] - xlim[0])*0.6, ylim[1] - (ylim[1] - ylim[0])*0.1,
             r'$1 / \Re(\lambda_2) = %.1f$ (years)' \
             % (-1. / eigvalGen[1].real,), fontsize=fs_latex)
    fig.savefig('%s/spectrum/eigval/figs/eigval_nev%d%s.%s' % (dstDir, nev, postfix, figFormat),
                bbox_inches='tight', dpi=dpi)
    
    # Plot eigenvectors of transfer operator
    tol = 0.
    alpha = 0.0
    for k in np.arange(nevPlot):
        if np.abs(eigval[k].real - 1) < 1.e-3:
            # Plot invariant measure
            print 'Plotting stationary density...'
            statDen = eigvec[0].real
            statDen /= statDen.sum()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            alpha = 0.01
            h = ax.contourf(X, Y, statDen.reshape(nx, nx), levels,
                            cmap=cm.hot_r)
            ax.set_xlim(gridXlim)
            ax.set_ylim(gridYlim)
            plt.colorbar(h)
            ax.set_xlabel("X Axis")
            ax.set_ylabel("Y Axis")
            ax.set_title("Approximation of the invariant measure", fontsize=fs_default)
            fig.savefig('%s/spectrum/eigvec/figs/eigvecReal_nev%d_ev03%d%s.%s' \
                        % (dstDir, nev, 1, postfix, figFormat), bbox_inches='tight', dpi=dpi)
        else:
            print 'Plotting real part of eigenvector %d...' % (k+1,)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            v2Real = eigvec[k].real
#            vmax = np.sort(np.abs(v2Real))[(1. - alpha)*N-1]
            vmax = np.sort(np.abs(v2Real))[-1]
            v2Real[v2Real > vmax] = vmax
            v2Real[v2Real < -vmax] = -vmax
            h = ax.contourf(X, Y, v2Real.reshape(nx, nx), levels,
                            cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
            ax.set_xlim(gridXlim)
            ax.set_ylim(gridYlim)
            plt.colorbar(h)
            ax.set_xlabel("X Axis")
            ax.set_ylabel("Y Axis")
            ax.set_title("Real part of the eigenvector %d" % (k+1,),
                         fontsize=fs_default)
            fig.savefig('%s/spectrum/eigvec/figs/eigvecReal_nev%d_ev03%d%s.%s' \
                        % (dstDir, nev, k+1, postfix, figFormat), bbox_inches='tight', dpi=dpi)

            if eigval[k].imag != 0:
                print 'Plotting imaginary  part of eigenvector %d...' % (k+1,)
                fig = plt.figure()
                ax = fig.add_subplot(111)
                v2Imag = eigvec[k].imag
                vmax = np.sort(np.abs(v2Imag))[-1]
                v2Imag[v2Imag > vmax] = vmax
                v2Imag[v2Imag < -vmax] = -vmax
                h = ax.contourf(X, Y, v2Imag.reshape(nx, nx), levels,
                                cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
                ax.set_xlim(gridXlim)
                ax.set_ylim(gridYlim)
                plt.colorbar(h)
                ax.set_xlabel("X Axis")
                ax.set_ylabel("Y Axis")
                ax.set_title("Imaginary part of the eigenvector %d" % (k+1,),
                             fontsize=fs_default)
                fig.savefig('%s/spectrum/eigvec/figs/eigvecImag_nev%d_ev03%d%s.%s' \
                            % (dstDir, nev, k, postfix, figFormat), bbox_inches='tight', dpi=dpi)

 
    # Plot eigenvectors of Koopman operator
    if plotAdjoint:
        eigvecAdjointScale = np.zeros((nevSingle, N), dtype=complex)
        for k in np.arange(nevPlot):
            eigvecAdjointScale[k] = eigvecAdjoint[k] \
                                    / np.conjugate(np.vdot(eigvecAdjoint[k], eigvec[k]))
            if np.abs(eigval[k].real - 1) < 1.e-3:
                # Plot invariant measure
                print 'Plotting ergodic vector...'
                ergodicVec = np.abs(eigvecAdjoint[0].real)
                ergodicVec /= np.abs(ergodicVec).max()
                fig = plt.figure()
                ax = fig.add_subplot(111)
                alpha = 0.01
                h = ax.contourf(X, Y, ergodicVec.reshape(nx, nx), levels,
                                cmap=cm.hot_r)
                ax.set_xlim(gridXlim)
                ax.set_ylim(gridYlim)
                plt.colorbar(h)
                ax.set_xlabel("X Axis")
                ax.set_ylabel("Y Axis")
                ax.set_title("Approximation of the ergodic vector", fontsize=fs_default)
                fig.savefig('%s/spectrum/eigvec/figs/eigvecAdjointReal_nev%d_ev03%d%s.%s' \
                            % (dstDir, nev, 1, postfix, figFormat), bbox_inches='tight', dpi=dpi)
            else:
                print 'Plotting real part of Koopman eigenvector %d...' % (k+1,)
                fig = plt.figure()
                ax = fig.add_subplot(111)
                v2Real = eigvecAdjoint[k].real
                vmax = np.sort(np.abs(v2Real))[-1]                        
                v2Real[v2Real > vmax] = vmax
                v2Real[v2Real < -vmax] = -vmax
                h = ax.contourf(X, Y, v2Real.reshape(nx, nx), levels,
                                cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
                ax.set_xlim(gridXlim)
                ax.set_ylim(gridYlim)
                plt.colorbar(h)
                ax.set_xlabel("X Axis")
                ax.set_ylabel("Y Axis")
                ax.set_title("Real part of the Koopman eigenvector %d" % (k+1,),
                             fontsize=fs_default)
                fig.savefig('%s/spectrum/eigvec/figs/eigvecAdjointReal_nev%d_ev03%d%s.%s' \
                            % (dstDir, nev, k+1, postfix, figFormat), bbox_inches='tight', dpi=dpi)

                if eigval[k].imag != 0:
                    print 'Plotting imaginary  part of Koopman eigenvector %d...' % (k+1,)
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    v2Imag = eigvecAdjoint[k].imag
                    vmax = np.sort(np.abs(v2Imag))[-1]
                    v2Imag[v2Imag > vmax] = vmax
                    v2Imag[v2Imag < -vmax] = -vmax
                    h = ax.contourf(X, Y, v2Imag.reshape(nx, nx), levels,
                                    cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
                    ax.set_xlim(gridXlim)
                    ax.set_ylim(gridYlim)
                    plt.colorbar(h)
                    ax.set_xlabel("X Axis")
                    ax.set_ylabel("Y Axis")
                    ax.set_title("Imaginary part of the Koopman eigenvector %d" % (k+1,),
                                 fontsize=fs_default)
                    fig.savefig('%s/spectrum/eigvec/figs/eigvecAdjointImag_nev%d_ev03%d%s.%s' \
                                % (dstDir, nev, k, postfix, figFormat),
                                bbox_inches='tight', dpi=dpi)

# Get ccf
lagMax = 100
lags = np.arange(0, lagMax+1, 1)
f = X.flatten()
g = X.flatten()
obsIdx0 = 0
obsIdx1 = 0

statesFileName = "%s/obs/obs%s.txt" % (dstDir, gridPostfix)
sim = np.loadtxt(statesFileName)
sim = sim.reshape(sim.shape[0] / dim, dim)
ccf = atmath.ccf(sim[:, obsIdx0], sim[:, obsIdx1], lagMax=lagMax)[lagMax:]
#ccf = atmath.ccovf(sim[:, obsIdx0], sim[:, obsIdx1], lagMax=lagMax)[lagMax:]

ccfRec = np.zeros((lags.shape[0],), dtype=complex)
for ev in np.arange(1, nevSingle):
    ccfRec += np.exp(eigvalGen[ev]*lags) \
              * (f * statDen * np.conjugate(eigvecAdjoint[ev])).sum() \
              * (eigvec[ev] * np.conjugate(g)).sum()
#ccfRec /= np.cov(sim[:, obsIdx0], sim[:, obsIdx1])[0, 1]
ccfRec /= ccfRec[0]

print 'Plotting'
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(lags, ccf, linewidth=2)
ax.plot(lags, ccfRec, '--', linewidth=2)
fig.savefig('%s/ccf/ccf_%s_%s_nev%d%s.%s' \
            % (dstDir, indexChoice[obsIdx0], indexChoice[obsIdx1], nev, postfix, figFormat),
            bbox_inches='tight', dpi=dpi)

    
