import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/helium_06/"
resDir = "../../../../res/plot/helium_06/";

fileNameBins1 = "bins_01.txt"
fileNameVariance1 = "energy_variance_01.txt"

fileNameBins2 = "bins_02.txt"
fileNameVariance2 = "energy_variance_02.txt"

fileNameBins3 = "bins_03.txt"
fileNameVariance3 = "energy_variance_03.txt"
# Load data.
os.chdir(curDir)
os.chdir(resDir)

blocks1  = np.loadtxt(fileNameBins1)
std1   = np.loadtxt(fileNameVariance1)

blocks2  = np.loadtxt(fileNameBins2)
std2   = np.loadtxt(fileNameVariance2)

blocks3  = np.loadtxt(fileNameBins3)
std3   = np.loadtxt(fileNameVariance3)

plt.plot(blocks1, std1)
plt.plot(blocks2, std2)
plt.plot(blocks3, std3)
plt.savefig('helium_06.png')
plt.show()
