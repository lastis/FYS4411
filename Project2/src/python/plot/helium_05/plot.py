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

font = {'family' : 'serif',
        'size'   : 15}

plt.rc('font',**font)

fig, ax = plt.subplots()

ax.plot(betaArray1, mean1,'b',label='GTO 3-21G')
ax.fill_between(betaArray1, mean1-meanVar1, mean1+meanVar1, color='b', alpha=0.2)


ax.set_ylabel(r'$E_0$ [a.u.]')
ax.set_xlabel(r'$\beta$ [a.u.]')
ax.legend(loc='best',fontsize=15)

ax.grid('on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('helium_05_pretty.png')
plt.show()


plt.plot(betaArray1, mean1)
plt.plot(betaArray1, mean1+meanVar1)
plt.plot(betaArray1, mean1-meanVar1)
plt.savefig('helium_05.png')
plt.show()
