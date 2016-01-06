import os
import numpy as np

SRng = np.array([1264, 1275, 1290, 1315, 1345, 1400])
restartState = 'warm'

#SRng = np.array([1400])
#restartState = "cold"
#restartState = "warm"

restartYear = 1100
lastYear = 9999
model = 'plasim'
rootDir = '../%s/' % model
exe = './most_plasim_salloc'
dataPrefix = 'MOST'
nproc = 1

cwd = os.getcwd()
for S in SRng:
    os.chdir(cwd)
    dstDir = '%s_%05d' % (restartState, S*10)
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
