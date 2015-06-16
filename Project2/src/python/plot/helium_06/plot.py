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



font = {'family' : 'serif',
        'size'   : 15}

plt.rc('font',**font)

fig, ax = plt.subplots()

ax.plot(blocks2, std2,'b',label=r'$\Delta t = 0.01$')
ax.plot(blocks1, std1,'r',label=r'$\Delta t = 0.001$')
ax.plot(blocks3, std3,'g',label=r'$\Delta t = 0.0001$')

ax.set_ylabel(r'Variance')
ax.set_xlabel(r'Block size')
ax.legend(loc='lower right', fontsize=15)

ax.grid('on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('helium_06_pretty.png')
plt.show()




# plt.plot(blocks1, std1)
# plt.plot(blocks2, std2)
# plt.plot(blocks3, std3)
# plt.savefig('helium_06.png')
# plt.show()
