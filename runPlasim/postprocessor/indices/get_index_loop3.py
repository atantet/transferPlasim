import os, sys
import numpy as np
from netCDF4 import Dataset

SRng = np.array([1267.5, 1272.5, 1285, 1295])
restartStateRng = ['warm']*SRng.shape[0]

#SRng = np.array([1272.5])
#restartStateRng = ['warm']

firstYear = 101
lastYear = 4000
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
