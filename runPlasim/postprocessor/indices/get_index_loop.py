import os, sys
import numpy as np
from netCDF4 import Dataset

#SRng = np.array([1260, 1360, 1380, 1400, 1415, 1425, 1430, 1433,
#                 1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
#restartStateRng = np.concatenate((['cold']*8, ['warm']*8), 0)

SRng = np.array([1263, 1265, 1270, 1280, 1300, 1330, 1360, 1435])
restartStateRng = ['warm']*8

#SRng = np.array([1435])
#restartStateRng = ['cold']
#restartStateRng = ['warm']

firstYear = 101
lastYear = 9999
yearsPerFile = 100
daysPerYear = 360

#indexChoice = 'globmst'
#indexChoice = 'npolemst'
#indexChoice = 'eqmst'
#indexChoice = 'eqlntr'
#indexChoice = 'energytransport'
#indexChoice = 'latmaxtransport'
indexChoice = 'globsic'
#indexChoice = 'nhemisic'
#indexChoice = 'globdep'
#indexChoice = 'eqdep'
#indexChoice = 'mtg'
script = 'get_index.py'

#indexChoice = 'areabelowtf20glob'
#indexChoice = 'areabelowtf20nhemi'
#indexChoice = 'areabelowtf10nhemi'
#indexChoice = 'areabelowtfnhemi'
#script = 'get_index_area.py'

for k in np.arange(SRng.shape[0]):
    S = SRng[k]
    restartState = restartStateRng[k]
    print 'Getting index %s for case %s_%d:' % (indexChoice, restartState, S)
    os.system('python %s %d %s %s %d %d %d %d' \
                  % (script, S, restartState, indexChoice,
                     firstYear, lastYear, yearsPerFile, daysPerYear))
