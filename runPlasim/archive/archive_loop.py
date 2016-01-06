import os
import numpy as np

#SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433])
#restartState = 'cold'

SRng = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartStateRng = ['warm']*8

#SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433,
#                 1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
#restartStateRng = np.concatenate((['cold']*8, ['warm']*8), 0)

prefixData = 'MOST'
win = 54   # win = 1 -> [101-200]
firstYear = 100 * win + 1
nyears = 100

archiveDir = '/archive/dijkbio/alexis/PlaSim/'
experimentDir = 'w2sb/'

os.system('mkdir %s/%s 2> /dev/null' % (archiveDir, experimentDir))
for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]

    srcDir = '../%s_%04d/' % (restartState, S)
    print 'Archiving %s from years %d to %d' \
        % (srcDir, firstYear, firstYear + nyears - 1)    

    # Make directory to be archived
    tarDir = '%s_%04d_%05d_%05d' \
        % (restartState, S, firstYear, firstYear + nyears - 1)
    os.system('mkdir %s/ 2> /dev/null' % tarDir)

    # Move files to directory to be archvied
    for year in np.arange(firstYear, firstYear + nyears):
        # Move data files
        recordFile = '%s.%03d' % (prefixData, year)
        os.system('mv %s/%s %s/' % (srcDir, recordFile, tarDir))
        # Move diagnostic files
        recordFile = '%s_DIAG.%03d' % (prefixData, year)
        os.system('mv %s/%s %s/' % (srcDir, recordFile, tarDir))
        # Move restart files but last year which is only copied
        recordFile = '%s_REST.%03d' % (prefixData, year)
        if year < firstYear + nyears - 1:
            os.system('mv %s/%s %s/' % (srcDir, recordFile, tarDir))
        else:
            os.system('cp %s/%s %s/' % (srcDir, recordFile, tarDir))

    # Create the tar-ball
    print 'Creating tar-ball %s.tar' % tarDir
    os.system('tar -cvf %s.tar %s' % (tarDir, tarDir))

    # Remove archived directory
    print 'Removing now-archived directory %s' % tarDir
    os.system('rm -r %s/' % tarDir)
    
    # Move to archive
    print 'Moving %s.tar to %s/%s/' % (tarDir, archiveDir, experimentDir)
    os.system('mv %s.tar %s/%s/' % (tarDir, archiveDir, experimentDir))

    print 'Done.'
