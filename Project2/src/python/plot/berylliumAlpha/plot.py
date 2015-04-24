import os
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
resDir = "../../../../res/plot/berylliumAlpha/";

os.chdir(resDir)
eneFile = open('energy.dat', 'r')
stdFile = open('energyStd.dat', 'r')
data1 = np.loadtxt(eneFile)
data2 = np.loadtxt(stdFile)
stdFile.close()
eneFile.close()

# plt.plot(data1)
# plt.show()

plt.plot(data2[0], data2[1],"x")
plt.show()

