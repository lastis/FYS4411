import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/berylliumSlaterQuantum/"
resDir = "../../../../res/plot/berylliumSlaterQuantum/";

alphaFileName = "alpha_array.txt"
fileName = "energies_mean.txt"

# Load data.
os.chdir(curDir)
os.chdir(resDir)
alphaFile   = open(alphaFileName,'r')
meanFile    = open(fileName, 'r')

alphaArray  = np.loadtxt(alphaFile)
meanArray   = np.loadtxt(meanFile)

alphaFile.close()
meanFile.close()

mean = np.mean(meanArray,axis=1)
meanVar = np.std(meanArray,axis=1)

plt.plot(alphaArray, mean)
plt.plot(alphaArray, mean+meanVar)
plt.plot(alphaArray, mean-meanVar)
plt.savefig('alpha_plot.png')
plt.show()
