import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/beryllium_02/"
resDir = "../../../../res/plot/beryllium_02/";

fileNameBins1 = "bins_01.txt"
fileNameVariance1 = "energy_variance_01.txt"

fileNameBins2 = "bins_02.txt"
fileNameVariance2 = "energy_variance_02.txt"

# Load data.
os.chdir(curDir)
os.chdir(resDir)

blocks1  = np.loadtxt(fileNameBins1)
std1   = np.loadtxt(fileNameVariance1)
# std1 = np.sqrt(std1)

blocks2  = np.loadtxt(fileNameBins2)
std2   = np.loadtxt(fileNameVariance2)
# std2 = np.sqrt(std2)


font = {'family' : 'serif',
        'size'   : 15}

plt.rc('font',**font)

fig, ax = plt.subplots()

ax.plot(blocks1, std1,'b',label='HL')

ax.set_ylabel(r'Variance')
ax.set_xlabel(r'Block size')
ax.legend(loc='lower right', fontsize=15)

ax.grid('on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('beryllium_02_pretty.png')
plt.show()



plt.plot(blocks1, std1)
# plt.plot(blocks2, std2)
plt.savefig('beryllium_02.png')
plt.show()
