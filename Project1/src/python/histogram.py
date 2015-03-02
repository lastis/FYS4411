import numpy as np
import matplotlib.pyplot as plt

infile = open('../../res/density.txt')

rMin = 0
rMax = int(infile.readline())

# density = infile.readlines()
density = np.loadtxt(infile)

l = len(density[0])
print l

p1 = density[0]

p2 = density[1]

pos = np.linspace(0,rMax,l)

plt.plot(pos, density[0])
plt.hold('on')
plt.plot(pos, density[1])
plt.show()
