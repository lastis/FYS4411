import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/helium_01/"
resDir = "../../../../res/plot/helium_01/";

betaFileName1 = "beta_array_01.txt"
fileName1 = "energies_01.txt";
betaFileName2 = "beta_array_02.txt"
fileName2 = "energies_02.txt";

# Load data.
os.chdir(curDir)
os.chdir(resDir)

betaArray1  = np.loadtxt(betaFileName1)
meanArray1   = np.loadtxt(fileName1)
mean1 = np.mean(meanArray1,axis=1)
meanVar1 = np.std(meanArray1,axis=1)

betaArray2  = np.loadtxt(betaFileName2)
meanArray2   = np.loadtxt(fileName2)
mean2 = np.mean(meanArray2,axis=1)
meanVar2 = np.std(meanArray2,axis=1)

print "Minimum 1: ", mean1.min() , " Value: ", betaArray1[mean1.argmin()]
print "Minimum 2: ", mean2.min() , " Value: ", betaArray2[mean2.argmin()]


font = {'family' : 'serif',
        'size'   : 15}

plt.rc('font',**font)

#hfont = {'fontname':'times'}
fig, ax = plt.subplots()

# ax.tick_params(axis='x', labelsize=18)
# ax.tick_params(axis='y', labelsize=18)
ax.plot(betaArray1, mean1,'r',label='With importance sampling')
ax.fill_between(betaArray1, mean1-meanVar1, mean1+meanVar1, color='r', alpha=0.2)
ax.plot(betaArray2, mean2,'b',label='Without importance sampling')
ax.fill_between(betaArray2, mean2-meanVar2, mean2+meanVar2, color='b', alpha=0.2)


ax.set_ylabel(r'$E_0$ [a.u.]')
ax.set_xlabel(r'$\beta$ [a.u.]')
ax.legend(loc='best',fontsize=15)

ax.grid('on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('helium_01_pretty.png')
plt.show()

# plt.plot(betaArray1, mean1)
# plt.plot(betaArray1, mean1+meanVar1)
# plt.plot(betaArray1, mean1-meanVar1)
# plt.plot(betaArray2, mean2)
# plt.plot(betaArray2, mean2+meanVar2)
# plt.plot(betaArray2, mean2-meanVar2)
# plt.savefig('helium_01.png')
# plt.show()
