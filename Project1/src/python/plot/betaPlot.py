from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from math import factorial

betafile = open('../../../res/heliumWave2Beta/beta.txt', 'r')
energyfile = open('../../../res/heliumWave2Beta/energy.txt', 'r')
variancefile = open('../../../res/heliumWave2Beta/variance.txt', 'r')
betaData = np.loadtxt(betafile)
energyData = np.loadtxt(energyfile)
varianceData = np.loadtxt(variancefile)
energyData = energyData.T
varianceData = varianceData.T

energyfile.close()
betafile.close()
variancefile.close()

# print len(energyData)
# print len(energyData[0])
energy = np.zeros(len(energyData))
energySq = np.zeros(len(energyData))
variance = np.zeros(len(energyData))
std = np.zeros(len(energyData))

N = len(energyData)
M = len(energyData[0])
print M

for i in xrange(N):
    for j in xrange(M):
	energy[i] += energyData[i][j]
        energySq[i] += energyData[i][j]*energyData[i][j]
    energy[i] = energy[i]/len(energyData[0])
    energySq[i] = energySq[i]/(len(energyData[0]))
    std[i] = np.sqrt(energySq[i] - energy[i]*energy[i])

plt.plot(betaData, energy)
plt.show()
plt.plot(betaData, std)
plt.show()

