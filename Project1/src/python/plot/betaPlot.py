from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

betafile = open('../../../res/heliumWave1Beta/beta.txt', 'r')
energyfile = open('../../../res/heliumWave1Beta/energy.txt', 'r')
betaData = np.loadtxt(betafile)
energyData = np.loadtxt(energyfile)
energyData = energyData.T

energyfile.close()
betafile.close()

# print len(energyData)
# print len(energyData[0])
energy = np.zeros(len(energyData))

N = len(energyData)
M = len(energyData[0])

for i in xrange(N):
    for j in xrange(M):
	energy[i] += energyData[i][j]
    energy[i] = energy[i]/len(energyData)

plt.plot(betaData, energy)

plt.show()

