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

energyArray1   = np.loadtxt(fileName1)
mean1 = np.mean(energyArray1,axis=0)
meanVar1 = np.std(energyArray1,axis=0)

energyArray2   = np.loadtxt(fileName2)
mean2 = np.mean(energyArray2,axis=0)
meanVar2 = np.std(energyArray2,axis=0)

font = {'family' : 'serif',
        'size'   : 15}

plt.rc('font',**font)

fig, ax = plt.subplots()

n = len(mean1)
x = range(0,n)

ax.plot(x, mean1,'r',label='HL')
ax.fill_between(x, mean1-meanVar1, mean1+meanVar1, color='r', alpha=0.2)
ax.plot(x, mean2,'b',label='GTO 3-21G')
ax.fill_between(x, mean2-meanVar2, mean2+meanVar2, color='b', alpha=0.2)

ax.set_ylabel(r'$E_0$ [a.u.]')
ax.set_xlabel(r'$\beta$')
ax.legend(loc='best',fontsize=15)

plt.axis([0,n,-30,0])

ax.grid('on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('beryllium_04_pretty.png')
plt.show()

# plt.plot(mean1)
# plt.plot(mean1+meanVar1)
# plt.plot(mean1-meanVar1)

# plt.plot(mean2)
# plt.plot(mean2+meanVar2)
# plt.plot(mean2-meanVar2)

# plt.savefig('beryllium_04.png')
# plt.show()
