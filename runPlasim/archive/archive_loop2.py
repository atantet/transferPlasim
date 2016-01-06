import os
import numpy as np

SRng = np.array([1264, 1275, 1290, 1315, 1345, 1400])
restartStateRng = ['warm']*6

prefixData = 'MOST'
win = 13   # win = 1 -> [101-200]
firstYear = 100 * win + 1
nyears = 100

archiveDir = '/archive/dijkbio/alexis/PlaSim/'
experimentDir = 'w2sb/'

os.system('mkdir %s/%s 2> /dev/null' % (archiveDir, experimentDir))
for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]

    srcDir = '../%s_%05d/' % (restartState, S*10)
    print 'Archiving %s from years %d to %d' \
        % (srcDir, firstYear, firstYear + nyears - 1)    

    # Make directory to be archived
    tarDir = '%s_%05d_%05d_%05d' \
        % (restartState, S*10, firstYear, firstYear + nyears - 1)
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
