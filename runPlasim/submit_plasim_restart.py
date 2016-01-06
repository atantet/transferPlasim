import os
import numpy as np

#SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433])
#restartState = 'cold'

SRng = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartState = 'warm'

#SRng = np.array([1330])
#restartState = "cold"
#restartState = "warm"

restartYear = 5800
lastYear = 9999
model = 'plasim'
rootDir = '../%s/' % model
exe = './most_plasim_salloc'
dataPrefix = 'MOST'
nproc = 1

cwd = os.getcwd()
for S in SRng:
    os.chdir(cwd)
    dstDir = '%s_%04d' % (restartState, S)
    # Go to working directory
    print 'cd to %s' % dstDir
    os.chdir(dstDir)

    # Copy proper restart file
    print 'Copying restart file %s_REST.%03d %s_restart' \
        % (dataPrefix, restartYear, model)
    os.system('cp %s_REST.%03d %s_restart' % (dataPrefix, restartYear, model))

    # Send job
    print 'Send job %s nproc=%d' % (exe, nproc)
    os.system('nohup %s %d %d %d > out.nohup & echo $! > nohup.pid' \
                  % (exe, restartYear, lastYear, nproc))
