from subprocess import call
import os
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
resDir = "../../../res/";
cppDir = "../../cpp/plot/"
os.chdir(cppDir)
call(["make","cpp_file=berylliumAlpha.cpp"])
call(["./a.out"])

os.chdir(curDir)
os.chdir(resDir)
stdFile = open('plot/berylliumAlpha/energyStd.dat', 'r')
data = np.loadtxt(stdFile)
stdFile.close()

plt.plot(data[0], data[1])
plt.show()

