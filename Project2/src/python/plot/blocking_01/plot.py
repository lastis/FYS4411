import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/blocking_01/"
resDir = "../../../../res/plot/blocking_01/";

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
std1 = np.sqrt(std1)

blocks2  = np.loadtxt(fileNameBins2)
std2   = np.loadtxt(fileNameVariance2)
std2 = np.sqrt(std2)

blocks3  = np.loadtxt(fileNameBins3)
std3   = np.loadtxt(fileNameVariance3)
std3 = np.sqrt(std3)

plt.plot(blocks1, std1)
plt.plot(blocks2, std2)
plt.plot(blocks3, std3)
plt.savefig('blocking_01.png')
plt.show()
