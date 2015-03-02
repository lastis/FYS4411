import numpy as np
import matplotlib.pyplot as plt

infile = open('../../res/density.txt')

density = np.loadtxt(infile)

l1 = len(density[1])
p1 = denisty[0][0]

l2 = len(density[3])
p2 = denisty[2][0]


plt.hist(0,p2,l2, density[1])
plt.hold('on')
plt.hist(linspace(0,p2,l2), density[3])
plt.show()
