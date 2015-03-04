from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

alphafile = open('../../../res/heliumWave2Alpha/alpha_energy_variance.txt', 'r')
data = np.loadtxt(alphafile)

plt.plot(data[0], data[1])

plt.show()

alphafile.close()
