import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/blocking/neon_prep/"
resDir = "../../../../res/blocking/neon_prep/";

blockingFileName = "energies_blocking.txt"
blockingStdFileName = "energies_blocking_std.txt"

# Load data.
os.chdir(curDir)
os.chdir(resDir)
blocksFile   = open(blockingFileName,'r')
stdFile    = open(blockingStdFileName, 'r')

blocks  = np.loadtxt(blocksFile)
std   = np.loadtxt(stdFile)

blocksFile.close()
stdFile.close()

std = np.sqrt(std)

plt.plot(blocks, std)
plt.savefig('blocking.png')
plt.show()
