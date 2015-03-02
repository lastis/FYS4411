from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

betafile = open('../../res/betaPlot.txt', 'r')
data = np.loadtxt(betafile)

# fig = plt.figure()
plt.plot(data[0], data[1])

plt.show()

betafile.close()
