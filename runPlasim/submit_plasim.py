import os
import numpy as np

#SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433])
#restartState = 'cold'

SRng = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartState = 'warm'

restartYear = 100
lastYear = 9999
model = 'plasim'
rootDir = '../%s/' % model
exe = './most_plasim_salloc'
nproc = 1

cwd = os.getcwd()
for S in SRng:
    os.chdir(cwd)
    # Create new working directory
    dstDir = '%s_%04d' % (restartState, S)
    print 'Create new working directory %s' % dstDir
    os.system('rm -r %s 2> /dev/null' % dstDir)
    os.system('mkdir %s 2> /dev/null' % dstDir)

    # Copy initial files to directory
    print 'Copying %s/restart/init/ to %s/' % (rootDir, dstDir)
    os.system('cp -r %s/restart/init/* %s/' % (rootDir, dstDir))

    # Copy proper restart file
    print 'Copying restart file %s/restart/%s_restart_%04d_%s to %s/%s_restart' % (rootDir, model, restartYear, restartState, dstDir, model)
    os.system('cp %s/restart/%s_restart_%04d_%s %s/%s_restart' % (rootDir, model, restartYear, restartState, dstDir, model))

    # Go to working directory
    print 'cd to %s' % dstDir
    os.chdir(dstDir)

    # Change the solar constant
    print 'Changing solar constant to %.1f in planet_namelist' % S
#    os.system('cat planet_namelist | sed "s/^.* GSOL0       =.*$/ GSOL0       = %.1f, /" > planet_namelist' % S)
    os.system('echo " GSOL0       =  %.1f" >> planet_namelist' % S)
    os.system('echo " /END" >> planet_namelist')
                                    
    # Send job
    print 'Send job %s nproc=%d' % (exe, nproc)
    os.system('nohup %s %d %d %d > out.nohup & echo $! > nohup.pid' \
                  % (exe, restartYear, lastYear, nproc))
