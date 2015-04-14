from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

alphafile = open('../../../res/berylliumAlphaBeta/alpha.txt', 'r')
betafile = open('../../../res/berylliumAlphaBeta/beta.txt', 'r')
energyfile = open('../../../res/berylliumAlphaBeta/energy.txt', 'r')
alpha = np.loadtxt(alphafile)
beta = np.loadtxt(betafile)
energy = np.loadtxt(energyfile)

font = {'family' : 'serif',
        'size'   : 10}

plt.rc('font', **font)

fig2, ax2 = plt.subplots()

X,Y = np.meshgrid(alpha,beta)

cax = ax2.imshow(energy,extent=[alpha[0],alpha[-1],beta[0],beta[-1]],\
        vmin=-16,vmax=-10,interpolation='lanczos',cmap=cm.coolwarm)
cbar = fig2.colorbar(cax, ticks=[-16, -13, -10], orientation='horizontal')

# plt.imshow(energy,extent=[alpha[0],alpha[-1],beta[0],beta[-1]])

ax2.set_title(r'Ground State Energies of Beryllium as Function of $\alpha$ and $\beta$' '\n' 
    r'with Generic Energy Calculation')

ax2.set_xlabel(r'$\beta$', fontsize=14)
ax2.set_ylabel(r'$\alpha$', fontsize=14)

ax2.grid('on')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.savefig('../../../res/berylliumAlphaBeta/berylliumAlphaBeta.eps')
plt.show()
plt.figure()


# Z = energy

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# X, Y = np.meshgrid(alpha, beta)

# Z = np.transpose(Z)

# surf = ax.plot_surface(X, Y, Z, 
#         rstride=1, cstride=1, cmap=cm.coolwarm,
#         linewidth=0, antialiased=False)
# # ax.set_zlim(-1.01, 1.01)

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.show()

alphafile.close()
betafile.close()
energyfile.close()
