# import numpy as np
# import matplotlib.pyplot as plt

# # Define the observable
# srcDir = '../runPlasim/postprocessor/indices/'

# # SRng = [1265, 1264, 1263]
# # lastYearRng = [9999, 6600, 9999]
# SRng = [1265, 1270, 1275, 1280, 1285, 1290]
# lastYearRng = [9999, 9999, 9999, 9999, 4000, 9999]

# dim = len(SRng)
# restartStateRng = ['warm'] * dim
# varName = 'eqmst'
# #varName = 'nhemisic'
# #varName = 'energytransport'
# restartState = 'warm'
# firstYear = 101
# daysPerYear = 360
# sampling = daysPerYear

# simYearly = []
# for k in np.arange(dim):
#     caseName = '%s_%d' % (restartStateRng[k], SRng[k]*10)
#     period = '%05d_%05d' % (firstYear, lastYearRng[k])
#     fileName = '%s/%s/%s_%s_%s.txt' \
#                % (srcDir, caseName, varName, caseName, period)
#     print 'Reading %s' % fileName
#     sim = np.loadtxt(fileName)
#     simYearly.append(np.convolve(sim, np.ones((sampling,)) \
#                                  / sampling)[::sampling])
# time = np.arange(simYearly[0].shape[0])

plotL = 1000
fig = plt.figure()
ax = fig.add_subplot(111)
for k in np.arange(dim):
    ax.plot(time[:plotL], simYearly[k][:plotL], label='%d' % SRng[k])
#ax.set_ylim(289, 294)
ax.set_ylim(14, 15.2)
plt.legend(loc='bottom right')
fig.savefig('spinup_%s_%d.png' % (varName, lastYearRng[1]),
            bbox_inches='tight', dpi=300)
                


