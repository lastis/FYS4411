import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/"
resDir = "../../../../res/plot/berylliumAlpha/";

fileName = "energies_mean.txt"
alphaFileName = "alpha_array.txt"

# # Delete the exsiting output file. 
# os.chdir(curDir)
# os.chdir(resDir)
# os.remove(fileName)
# # Compile the code
# os.chdir(curDir)
# os.chdir(cppDir)
# call(["make","cpp_file=berylliumAlphas.cpp"])

# # Run the cpp code
# os.chdir(curDir)
# os.chdir(cppDir)
# for alpha in alphaArray:
#     call(["./a.out", fileName, str(alpha),"0.8", "100"])

mean = []
meanVar = []
# alphaArray = np.linspace(3.3,4,5)
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
    meanVar.append(np.sqrt(sampleSq - sample*sample))

mean = np.asarray(mean)
meanVar = np.asarray(meanVar)
plt.plot(alphaArray, mean)
plt.plot(alphaArray, mean+meanVar)
plt.plot(alphaArray, mean-meanVar)
plt.show()

