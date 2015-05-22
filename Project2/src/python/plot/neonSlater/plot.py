import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/neonSlater/"
resDir = "../../../../res/plot/neonSlater/";

alphaFileName = "alpha_array.txt"
fileName = "energies_mean.txt"
fileNameDD = "DD.txt"
fileNameCC = "CC.txt"
fileNameDC = "DC.txt"
fileNameEPot = "energies_mean_potential.txt"

# Load data.
os.chdir(curDir)
os.chdir(resDir)
alphaFile   = open(alphaFileName,'r')
meanFile    = open(fileName, 'r')
fileDD      = open(fileNameDD,'r')
fileCC      = open(fileNameCC,'r')
fileDC      = open(fileNameDC,'r')
fileEPot      = open(fileNameEPot,'r')

alphaArray  = np.loadtxt(alphaFile)
meanArray   = np.loadtxt(meanFile)
DD          = np.loadtxt(fileDD)
CC          = np.loadtxt(fileCC)
DC          = np.loadtxt(fileDC)
EPot          = np.loadtxt(fileEPot)

alphaFile.close()
meanFile.close()
fileDD.close();
fileCC.close();
fileDC.close();
fileEPot.close();

mean = np.mean(meanArray,axis=1)
meanVar = np.std(meanArray,axis=1)
meanDD = np.mean(DD,axis=1)
meanCC = np.mean(CC,axis=1)
meanDC = np.mean(DC,axis=1)
meanEPot = np.mean(EPot,axis=1)

plt.plot(alphaArray, meanEPot);
plt.plot(alphaArray, meanDD);
plt.plot(alphaArray, meanCC);
plt.plot(alphaArray, meanDC);
plt.savefig('special.png')
plt.show()
plt.figure()

plt.plot(alphaArray, mean)
plt.plot(alphaArray, mean+meanVar)
plt.plot(alphaArray, mean-meanVar)
plt.savefig('alpha_plot.png')
plt.show()
