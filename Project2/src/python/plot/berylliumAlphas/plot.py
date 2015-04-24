
import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
cppDir = "../../../cpp/plot/"
resDir = "../../../../res/plot/berylliumAlpha/";

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=berylliumAlphas.cpp"])

mean = []
meanVar = []
alphaArray = np.linspace(3.3,4,100)
for alpha in alphaArray:
    # Run the cpp code
    os.chdir(curDir)
    os.chdir(cppDir)
    call(["./a.out", str(alpha),"0.8", "100"])

    # Manipulate data
    os.chdir(curDir)
    os.chdir(resDir)
    meanFile = open('array_mean.dat', 'r')
    meanArray = np.loadtxt(meanFile)
    meanFile.close()
    sample = 0
    sampleSq = 0
    for tmp in meanArray:
        sample += tmp
        sampleSq += tmp*tmp
    sample = sample/len(meanArray)
    sampleSq = sampleSq/len(meanArray)
    mean.append(sample)
    meanVar.append(sampleSq - sample*sample)

mean = np.asarray(mean)
meanVar = np.asarray(meanVar)
plt.plot(alphaArray, mean)
plt.plot(alphaArray, mean+meanVar)
plt.plot(alphaArray, mean-meanVar)
plt.show()

