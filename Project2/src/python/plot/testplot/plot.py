import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/testplot/"
resDir = "../../../../res/plot/testplot/";

fileName = "energies_mean.txt"
alphaFileName = "alpha_array.txt"

fileNameDD = "DD.txt"
fileNameCC = "CC.txt"
fileNameDC = "DC.txt"

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

fileDD = open(fileNameDD,'r')
fileCC = open(fileNameCC,'r')
fileDC = open(fileNameDC,'r')
DD = np.loadtxt(fileDD)
CC = np.loadtxt(fileCC)
DC = np.loadtxt(fileDC)
meanDD = []
meanCC = []
meanDC = []


if len(meanArray.shape) == 1:
    plt.plot(alphaArray, meanArray)
    plt.savefig('alpha_plot.png')
    plt.show()
else :
    for i,array in enumerate(DD):
        sample = 0
        for tmp in DD[i]:
            sample += tmp
        sample = sample/len(DD[i])
        meanDD.append(sample)
    meanDD = np.asarray(meanDD)

    for i,array in enumerate(CC):
        sample = 0
        for tmp in CC[i]:
            sample += tmp
        sample = sample/len(CC[i])
        meanCC.append(sample)
    meanCC = np.asarray(meanCC)

    for i,array in enumerate(DC):
        sample = 0
        for tmp in DC[i]:
            sample += tmp
        sample = sample/len(DC[i])
        meanDC.append(sample)
    meanDC = np.asarray(meanDC)

    for i,array in enumerate(meanArray):
        sample = 0
        sampleSq = 0
        for tmp in meanArray[i]:
            sample += tmp
            sampleSq += tmp*tmp

        sample = sample/len(meanArray[i])
        sampleSq = sampleSq/len(meanArray[i])
        mean.append(sample)
        # Actually the standard deviance of the mean.
        meanVar.append(np.sqrt((sampleSq - sample*sample)/len(meanArray[i])))
    mean = np.asarray(mean)
    meanVar = np.asarray(meanVar)



    plt.plot(alphaArray, meanDD);
    plt.plot(alphaArray, meanCC);
    plt.plot(alphaArray, meanDC);
    plt.savefig('special.png')
    plt.show()
    plt.figure()

    plt.plot(alphaArray, mean)
    plt.plot(alphaArray, mean+meanVar)
    plt.plot(alphaArray, mean-meanVar)
    plt.savefig('alpha_plot.png')
    plt.show()
