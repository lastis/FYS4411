import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/beryllium_05/"
resDir = "../../../../res/plot/beryllium_05/";

fileName1 = "density_01.txt";

# Load data.
os.chdir(curDir)
os.chdir(resDir)

densityMat1   = np.loadtxt(fileName1)
densityMean1 = np.mean(densityMat1,axis=0)
densityStd1 = np.std(densityMat1,axis=0)

rMax = 5
bins1 = np.linspace(0,rMax,len(densityMean1))

plt.plot(bins1, densityMean1)
plt.plot(bins1, densityMean1+densityStd1)
plt.plot(bins1, densityMean1-densityStd1)
plt.savefig('beryllium_05.png')
plt.show()
