import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/beryllium_03/"
resDir = "../../../../res/plot/beryllium_03/";

fileName1 = "density_01.txt";
fileName2 = "density_02.txt";

# Load data.
os.chdir(curDir)
os.chdir(resDir)

densityMat1   = np.loadtxt(fileName1)
densityMean1 = np.mean(densityMat1,axis=0)
densityStd1 = np.std(densityMat1,axis=0)

densityMat2   = np.loadtxt(fileName2)
densityMean2 = np.mean(densityMat2,axis=0)
densityStd2 = np.std(densityMat2,axis=0)

rMax = 5
bins1 = np.linspace(0,rMax,len(densityMean1))
bins2 = np.linspace(0,rMax,len(densityMean2))

plt.plot(bins1, densityMean1)
plt.plot(bins1, densityMean1+densityStd1)
plt.plot(bins1, densityMean1-densityStd1)
plt.plot(bins2, densityMean2)
plt.plot(bins2, densityMean2+densityStd2)
plt.plot(bins2, densityMean2-densityStd2)
plt.savefig('beryllium_03.png')
plt.show()
