import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/berylliumAlpha/"
resDir = "../../../../res/plot/berylliumAlpha/";

fileName = "energies_mean.txt"
alphaFileName = "alpha_array.txt"

mean = []
meanVar = []
# Manipulate data
os.chdir(curDir)
os.chdir(resDir)

alphaFile = open(alphaFileName,'r')
alphaArray = np.loadtxt(alphaFile)
alphaFile.close()

meanFile = open(fileName, 'r')
meanArray = np.loadtxt(meanFile)
meanFile.close()
for i,array in enumerate(meanArray):
    sample = 0
    sampleSq = 0
    for tmp in meanArray[i]:
        sample += tmp
        sampleSq += tmp*tmp

    sample = sample/len(meanArray[i])
    sampleSq = sampleSq/len(meanArray[i])
    mean.append(sample)
    # meanVar.append(np.sqrt(sampleSq - sample*sample))
    # Actually the SEM, Standard Error of the Mean. Or is it?
    meanVar.append(np.sqrt((sampleSq - sample*sample)/len(meanArray[i])))

mean = np.asarray(mean)
meanVar = np.asarray(meanVar)
plt.plot(alphaArray, mean)
plt.plot(alphaArray, mean+meanVar)
plt.plot(alphaArray, mean-meanVar)
plt.savefig('alpha_plot.png')
plt.show()
