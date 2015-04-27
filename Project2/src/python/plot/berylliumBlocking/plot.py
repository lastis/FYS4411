import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/berylliumBlocking/"
resDir = "../../../../res/plot/berylliumBlocking/";

fileName = "std_vs_blocksize.txt"

os.chdir(curDir)
os.chdir(resDir)

dataFile = open(fileName,'r')
data = np.loadtxt(dataFile)
dataFile.close()

plt.plot(data[1], data[0], 'x')
plt.show()
