from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

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
energySpline = np.zeros(len(energyData))
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

spline = interpolate.splrep(betaData,energy, s=0.00006)
betaNew = np.linspace(betaData[0],betaData[-1],1001)
energyNew = interpolate.splev(betaNew,spline)

minIndex = np.argmin(energyNew)
minBeta = betaNew[minIndex]
print minBeta


font = {'family' : 'serif',
        'size'   : 10}

plt.rc('font', **font)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(betaData, energy,'b')
ax.fill_between(betaNew, energyNew-std[0], energyNew+std[0], color='c', alpha=0.2)

ax.set_title(r'Ground State Energies as Function of $\beta$ for Helium with' '\n' 
    r'Local Energy Calculation and Trial Wavefunction $\psi_{T2}$')
ax.set_xlabel(r'$\beta$', fontsize=14)
ax.set_ylabel(r'$E_0$', fontsize=14)
ax.set_xlim(0,2.5)
ax.set_ylim(-2.88,-2.78)

ax.grid('on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('../../../res/heliumWave2Beta/heliumWave2Beta.eps')
plt.show()
