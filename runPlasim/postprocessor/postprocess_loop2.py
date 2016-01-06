import os
import numpy as np

#SRng = np.array([1264, 1275, 1290, 1315, 1345, 1400])
#restartStateRng = ['warm']*6
#variables = ['tsa', 'sic', 'ntr', 'gttd', 'hfns']
#win = 13

SRng = np.array([1345])
restartStateRng = ['warm']
variables = ['tsa']
win = 10

postprocessor = 'burn7dotYYYY.x'
sbatchScript = 'submit_postproc.sh'
prefixData = 'MOST.'
firstYear = 100 * win + 1
nyears = 100

dstDir = 'w2sb'
# Make post-processing directory
os.system('mkdir %s 2> /dev/null' % dstDir)
for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]
    srcDir = '../%s_%05d/' % (restartState, S*10)
    print 'Postprocessing %s from %d to %d' % (srcDir, firstYear, firstYear + nyears - 1)

    # Send job
    firstFile = '%s%03d' % (prefixData, firstYear)
    for ivar in np.arange(len(variables)):
        namelistFile = 'burn_namelist_%s%d' % (variables[ivar], nyears)
        dstFile = '%s_%s_%05d_%05d_%05d.nc' \
            % (variables[ivar], restartState, S*10,
               firstYear, firstYear + nyears - 1)
        print 'Burning with namelist %s' % namelistFile
        os.system('sbatch %s %s %s/%s %s/%s %s' \
                      % (sbatchScript, postprocessor, srcDir, firstFile,
                         dstDir, dstFile, namelistFile))
