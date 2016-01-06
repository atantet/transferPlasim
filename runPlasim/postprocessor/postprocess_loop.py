import os
import numpy as np

#SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433])
#restartStateRng = ['cold']*8
SRng = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartStateRng = ['warm']*8
#SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433,
#                 1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
#restartStateRng = np.concatenate((['cold']*8, ['warm']*8), 0)
variables = ['tsa', 'sic', 'ntr', 'gttd', 'hfns']

#SRng = np.array([1415])
#restartStateRng = ['cold']
#restartStateRng = ['warm']
#variables = ['gttd']

win = 60

postprocessor = 'burn7dotYYYY.x'
sbatchScript = 'submit_postproc.sh'
prefixData = 'MOST.'
   # win = 1 -> [101-200]
firstYear = 100 * win + 1
nyears = 100

dstDir = 'w2sb'
# Make post-processing directory
os.system('mkdir %s 2> /dev/null' % dstDir)
for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]
    srcDir = '../%s_%04d/' % (restartState, S)
    print 'Postprocessing %s from %d to %d' % (srcDir, firstYear, firstYear + nyears - 1)

    # Send job
    firstFile = '%s%03d' % (prefixData, firstYear)
    for ivar in np.arange(len(variables)):
        namelistFile = 'burn_namelist_%s%d' % (variables[ivar], nyears)
        dstFile = '%s_%s_%04d_%05d_%05d.nc' \
            % (variables[ivar], restartState, S,
               firstYear, firstYear + nyears - 1)
        print 'Burning with namelist %s' % namelistFile
        os.system('sbatch %s %s %s/%s %s/%s %s' \
                      % (sbatchScript, postprocessor, srcDir, firstFile,
                         dstDir, dstFile, namelistFile))
