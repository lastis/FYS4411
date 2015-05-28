import os
import numpy as np
# from mayavi.mlab import *
from numpy.random import uniform, seed
import scipy.interpolate
import matplotlib.pyplot as plt

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/visualization/neon/"

resDir = "../../../../res/visualization/neon/";
fileName = "density3D.txt"
# Read output from file
os.chdir(curDir)
os.chdir(resDir)

x,y,z,f  = np.loadtxt(fileName,unpack=True)

f = f/f.max()
points = np.asarray((x,y))
points = points.T
values = f

N = 101;
xMin = x.min();
yMin = y.min();
zMin = z.min();
xMax = np.array([x.max(),abs(xMin)])
xMax = xMax.max()
yMax = np.array([y.max(),abs(yMin)])
yMax = yMax.max()
zMax = np.array([z.max(),abs(zMin)])
zMax = zMax.max()

rMax = np.sqrt(xMax**2 + yMax**2 + zMax**2)
xArr = np.linspace(xMin,xMax,N)
yArr = np.linspace(yMin,yMax,N)
zArr = np.linspace(zMin,zMax,N)


xi, yi, zi = np.linspace(x.min(), x.max(), N), np.linspace(y.min(), y.max(), N), np.linspace(z.min(), z.max(), N)

xGrid, yGrid = np.meshgrid(xi, yi)

X = np.zeros([N,N])
Y = np.zeros([N,N])
Z = np.zeros([N,N])
density = np.zeros(N+1)

cnt=0
for i in range(0,N):
    for j in range(0,N):
        for k in range(0,N):
            X[j,k] += f[cnt]
            Y[i,k] += f[cnt]
            Z[i,j] += f[cnt]
            rAbs = np.sqrt(xArr[i]**2+yArr[j]**2+zArr[k]**2)
            rBin = int(rAbs*N/rMax)
            density[rBin] += f[cnt]
            cnt += 1


Z = scipy.interpolate.griddata(points,values,(xGrid,yGrid), method='nearest')

r = np.zeros(len(x))
plt.plot(density)
for i in range(len(x)):
    r[i] = np.sqrt(x[i]**2+y[i]**2+z[i]**2)

fig, axs = plt.subplots(3, sharex=True, sharey=True)
cs = axs[0].contourf(yi, zi, X)
fig.colorbar(cs, ax=axs[0], format="%.2f")

cs = axs[1].contourf(xi, zi, Y)
fig.colorbar(cs, ax=axs[1], format="%.2f")

cs = axs[2].contourf(xi, yi, Z)
fig.colorbar(cs, ax=axs[2], format="%.2f")
plt.show()

plt.imshow(X,interpolation = 'nearest')
plt.colorbar(format="%.2f")
plt.show()
