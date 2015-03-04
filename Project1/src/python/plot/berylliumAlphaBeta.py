from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

alphafile = open('../../../res/berylliumAlphaBeta/alpha.txt', 'r')
betafile = open('../../../res/berylliumAlphaBeta/beta.txt', 'r')
energyfile = open('../../../res/berylliumAlphaBeta/energy.txt', 'r')
alpha = np.loadtxt(alphafile)
beta = np.loadtxt(betafile)
energy = np.loadtxt(energyfile)
Z = energy

fig = plt.figure()
ax = fig.gca(projection='3d')

X, Y = np.meshgrid(alpha, beta)

Z = np.transpose(Z)

surf = ax.plot_surface(X, Y, Z, 
        rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

alphafile.close()
betafile.close()
energyfile.close()
