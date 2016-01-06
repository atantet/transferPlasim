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
indexChoice = ["nhemisic"]
dim = len(indexChoice)

# Grid definition
nx = 50
nSTD = 5

# Lags
nLags = 1
tauDimRng = np.array([1])

nev = 10
nevPlot = 1
#plotAdjoint = False
plotAdjoint = True
getCondition = True
plotCCF = True

resDir = "%s_%d" % (restartState, np.round(S*10))
dstDir = resDir
caseName = "%s_%05d_%05d" % (resDir, firstYear, lastYear)
os.system('mkdir %s/ccf %s/spectrum/condition %s/spectrum/eigval/figs %s/spectrum/eigvec/figs 2> /dev/null' % (dstDir, dstDir, dstDir, dstDir))

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
x = coord[0]

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
gridXlim = [coord[0].min(), coord[0].max()]


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

    if plotAdjoint or getCondition or plotCCF:
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

    # Get generator eigenvalues
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

    # Get condition number
    print 'Getting condition number...'
    conditionEigval = np.empty((nevSingle,))
    if getCondition:
#        normEigvec = np.linalg.cond(eigvec.T, p=2)
#        normEigvecAdjoint = np.linalg.cond(eigvecAdjoing.T, p=2)
#        conditionNumber = normEigvec * normEigvecAdjoint
 #       np.savetxt('%s/spectrum/condition/conditionNumber_nev%d%s.txt' \
 #                  % (dstDir, nev, postfix), (conditionNumber,))
        for ev in np.arange(nevSingle):
            conditionEigval[ev] = np.abs(np.sum(eigvec[ev] * np.conjugate(eigvec[ev]))) * np.abs(np.sum(eigvecAdjoint[ev] * np.conjugate(eigvecAdjoint[ev]))) / np.abs(np.sum(eigvec[ev] * np.conjugate(eigvecAdjoint[ev])))
        np.savetxt('%s/spectrum/condition/conditionEigval_nev%d%s.txt' \
                   % (dstDir, nev, postfix), conditionEigval)
    
    # Plot spectrum
    print 'Plotting spectrum slowest rate ', -1. / eigvalGen[1].real
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(eigvalGen[1:].real, eigvalGen[1:].imag,
               c='k', s=msize, marker='o')
    ax.scatter(eigvalGen[0].real, eigvalGen[0].imag,
               c='r', s=msize, marker='o')
    ax.set_xlabel(r'$\Re(\zeta_i)$', fontsize=fs_latex)
    ax.set_ylabel(r'$\Im(\zeta_i)$', fontsize=fs_latex)
    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    #ax.set_title('%d-time-step spectrum for %s\nSlowest time-scale: %.1f' \
        #    % (tau, srcPostfix, -1. / rate[0]))
    ax.set_xlim(-1.2, 0.02)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.text(xlim[1] - (xlim[1] - xlim[0])*0.6,
             ylim[1] - (ylim[1] - ylim[0])*0.1,
             r'$1 / \Re(\lambda_2) = %.1f$ (years)' \
             % (-1. / eigvalGen[1].real,), fontsize=fs_latex)
    fig.savefig('%s/spectrum/eigval/figs/eigval_nev%d%s.%s' \
                % (dstDir, nev, postfix, figFormat),
                bbox_inches='tight', dpi=dpi)

    # Plot condition number for each eigenvalue
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(1, nevSingle + 1)[:5], conditionEigval[:5], '+k')
#    ax.set_title('Condition number = %f' % conditionNumber)
    fig.savefig('%s/spectrum/condition/conditionEigval_nev%d%s.%s' \
                % (dstDir, nev, postfix, figFormat),
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
            h = ax.fill(x, statDen)
            ax.set_xlim(gridXlim)
            ax.set_xlabel("X Axis")
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
            h = ax.fill(x, v2Real)
            ax.set_xlim(gridXlim)
            ax.set_xlabel("X Axis")
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
                h = ax.fill(x, v2Imag)

                ax.set_xlim(gridXlim)
                ax.set_xlabel("X Axis")
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
                h = ax.fill(x, ergodicVec)
                ax.set_xlim(gridXlim)
                ax.set_xlabel("X Axis")
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
                h = ax.fill(x, v2Real)
                ax.set_xlim(gridXlim)
                ax.set_xlabel("X Axis")
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
                    h = ax.fill(x, v2Imag)
                    ax.set_xlim(gridXlim)
                    ax.set_xlabel("X Axis")
                    ax.set_title("Imaginary part of the Koopman eigenvector %d" % (k+1,),
                                 fontsize=fs_default)
                    fig.savefig('%s/spectrum/eigvec/figs/eigvecAdjointImag_nev%d_ev03%d%s.%s' \
                                % (dstDir, nev, k, postfix, figFormat),
                                bbox_inches='tight', dpi=dpi)

if plotCCF:
    # Get ccf
    lagMax = 100
    lags = np.arange(0, lagMax+1, 1)
    f = x
    g = f
    obsIdx0 = 0
    obsIdx1 = 0
    
    acfRec = np.zeros((lags.shape[0],), dtype=complex)
    for ev in np.arange(1, nevSingle):
        acfRec += np.exp(eigvalGen[ev]*lags) \
                  * (f * statDen * np.conjugate(eigvecAdjoint[ev])).sum() \
                  * (eigvec[ev] * np.conjugate(g)).sum()
        print (f * statDen * np.conjugate(eigvecAdjoint[ev])).sum() * (eigvec[ev] * np.conjugate(g)).sum()
    acfRec /= acfRec[0]

    statesFileName = "%s/obs/obs%s.txt" % (dstDir, gridPostfix)
    sim = np.loadtxt(statesFileName)
    sim = sim.reshape(sim.shape[0] / dim, dim)
    acf = atmath.ccf(sim[:, obsIdx0], sim[:, obsIdx1], lagMax=lagMax)[lagMax:]
    
    print 'Plotting'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(lags, acf, linewidth=2)
    ax.plot(lags, acfRec, '--', linewidth=2)
    ax.set_xlabel('lag (yr)', fontsize=fs_default)
                  ax.set_ylabel(r'$\Re(\lambda_i)$ (%s$^{-1}$)' % unitTime, fontsize=fs_latex)
    ax.set_xlim(xlim)
    ax.set_xticks(SRng)
    ax.set_xticklabels(xticklabels)
    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    ax.set_ylabel(r'$\Re(\lambda_i)$ (%s$^{-1}$)' % unitTime, fontsize=fs_latex)
    ax.set_ylim(ylimCCF)
    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    fig.savefig('%s/ccf/ccf_%s_%s_nev%d%s.%s' \
                % (dstDir, indexChoice[obsIdx0], indexChoice[obsIdx1],
                   nev, postfix, figFormat),
                bbox_inches='tight', dpi=dpi)

    
