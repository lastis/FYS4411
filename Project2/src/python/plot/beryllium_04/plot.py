import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/beryllium_04/"
resDir = "../../../../res/plot/beryllium_04/";

fileName1 = "energies_01.txt";
fileName2 = "energies_02.txt";

# Load data.
os.chdir(curDir)
os.chdir(resDir)

# energyArray1   = np.loadtxt(fileName1)
# mean1 = np.mean(energyArray1,axis=0)
# meanVar1 = np.std(energyArray1,axis=0)

energyArray2   = np.loadtxt(fileName2)
mean2 = np.mean(energyArray2,axis=0)
meanVar2 = np.std(energyArray2,axis=0)

# plt.plot(mean1)
# plt.plot(mean1+meanVar1)
# plt.plot(mean1-meanVar1)

plt.plot(mean2)
plt.plot(mean2+meanVar2)
plt.plot(mean2-meanVar2)

plt.savefig('beryllium_04.png')
plt.show()
