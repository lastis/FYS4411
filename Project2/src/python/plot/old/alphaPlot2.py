from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

alphafile = open('../../../res/heliumWave2Alpha/alpha_energy_variance.txt', 'r')
data = np.loadtxt(alphafile)

minIndex = np.argmin(data[1])
minAlpha = data[0][minIndex]
print minAlpha
plt.plot(data[0], data[1])
plt.ylim([0,2])

plt.show()

alphafile.close()
