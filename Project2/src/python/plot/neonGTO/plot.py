import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/neonGTO/"
resDir = "../../../../res/plot/neonGTO/";

betaFileName = "beta_array.txt"
fileName = "energies_mean.txt"

# Load data.
os.chdir(curDir)
os.chdir(resDir)
betaFile   = open(betaFileName,'r')
meanFile    = open(fileName, 'r')

betaArray  = np.loadtxt(betaFile)
meanArray   = np.loadtxt(meanFile)

betaFile.close()
meanFile.close()

mean = np.mean(meanArray,axis=1)
meanVar = np.std(meanArray,axis=1)

plt.plot(betaArray, mean)
plt.plot(betaArray, mean+meanVar)
plt.plot(betaArray, mean-meanVar)
plt.savefig('beta_plot.png')
plt.show()
