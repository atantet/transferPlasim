import os, sys
import numpy as np
from netCDF4 import Dataset

#SRng = np.array([1264, 1275, 1290, 1315, 1345, 1400])
#restartStateRng = ['warm']*6

SRng = np.array([1262.8, 1262.9, 1263.5, 1264.5])
restartStateRng = ['warm']*SRng.shape[0]

#SRng = np.array([1262.8])
#restartStateRng = ['warm']

firstYear = 101
lastYear = 6300
yearsPerFile = 100
daysPerYear = 360

#indexChoice = 'globmst'
#indexChoice = 'npolemst'
#indexChoice = 'eqmst'
#indexChoice = 'mtg'
#indexChoice = 'npolentr'
#indexChoice = 'energytransport'
#indexChoice = 'latmaxtransport'
#indexChoice = 'globsic'
#indexChoice = 'nhemisic'
#indexChoice = 'globdep'
indexChoice = 'eqdep'
script = 'get_index.py'

for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]
    print 'Getting index %s for case %s_%d:' \
        % (indexChoice, restartState, S*10)
    os.system('python %s %d %s %s %d %d %d %d' \
                  % (script, S*10, restartState, indexChoice,
                     firstYear, lastYear, yearsPerFile, daysPerYear))
