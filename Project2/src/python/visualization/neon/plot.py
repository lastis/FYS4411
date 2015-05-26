import os
import numpy as np
from mayavi.mlab import *
from numpy.random import uniform, seed
from scipy.interpolate import griddata

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/visualization/neon/"

resDir = "../../../../res/visualization/neon/";
fileName = "density3D.txt"
# Read output from file
os.chdir(curDir)
os.chdir(resDir)

x,y,z,f  = np.loadtxt(fileName,unpack=True)
points = np.asarray((x,y,z))
points = points.T
values = f


N = 201;
xMin = x[0];
yMin = y[0];
zMin = z[0];
xMax = x[-1];
yMax = y[-1];
zMax = z[-1];
xArr = np.linspace(xMin,xMax,N)
yArr = np.linspace(yMin,yMax,N)
zArr = np.linspace(zMin,zMax,N)

xGrid, yGrid, zGrid = np.meshgrid(xArr,yArr,zArr)

F = griddata(points,values,(xGrid,yGrid, zGrid), method='nearest')

contour3d(F,contours=20,opacity=.2)
show()
