import os
import numpy as np
import matplotlib.pyplot as plt

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/visualization/neon/"

resDir = "../../../../res/visualization/neon/";
fileName = "density1D.txt"
# Read output from file
os.chdir(curDir)
os.chdir(resDir)

r,density  = np.loadtxt(fileName,unpack=True)


N = 100;
rMin = r[0];
rMax = r[-1];
r = r*0.53

plt.plot(r,density)
plt.savefig('density1d.png')
plt.show()
