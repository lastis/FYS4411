import numpy as np
import matplotlib.pyplot as plt

infile = open('../../res/density.txt')

density = np.loadtxt(infile)

plt.hist(density[0])
plt.show()
plt.hist(density[1])
plt.show()
