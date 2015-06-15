from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/beryllium_05/"
resDir = "../../../../res/plot/beryllium_05/";

# Load data.
os.chdir(curDir)
os.chdir(resDir)

alphaFileName1 = "alpha_array_01.txt"
alphaArray1  = np.loadtxt(alphaFileName1)

betaFileName1 = "beta_array_01.txt"
betaArray1  = np.loadtxt(betaFileName1)

energyFileName1 = "energy_array_01.txt";
energyArray1  = np.loadtxt(energyFileName1)

X, Y = np.meshgrid(alphaArray1, betaArray1)
Z = energyArray1.T

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=cm.coolwarm)
# ax.set_zlim(-3,0)


plt.savefig('beryllium_05.png')
plt.show()
