import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/helium_05/"
resDir = "../../../../res/plot/helium_05/";

betaFileName1 = "beta_array_01.txt"
fileName1 = "energies_01.txt";

# Load data.
os.chdir(curDir)
os.chdir(resDir)

betaArray1  = np.loadtxt(betaFileName1)
meanArray1   = np.loadtxt(fileName1)
mean1 = np.mean(meanArray1,axis=1)
meanVar1 = np.std(meanArray1,axis=1)

print "Minimum 1: ", mean1.min() , " Value: ", betaArray1[mean1.argmin()]

plt.plot(betaArray1, mean1)
plt.plot(betaArray1, mean1+meanVar1)
plt.plot(betaArray1, mean1-meanVar1)
plt.savefig('helium_05.png')
plt.show()
