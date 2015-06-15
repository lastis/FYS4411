import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/ground_state/neon_01/"
resDir = "../../../../res/ground_state/neon_01/";

fileName1 = "energies_01.txt";
# fileName2 = "energies_02.txt";

# Load data.
os.chdir(curDir)
os.chdir(resDir)

meanArray1   = np.loadtxt(fileName1)
mean1 = np.mean(meanArray1)
meanVar1 = np.std(meanArray1)

# meanArray2   = np.loadtxt(fileName2)
# mean2 = np.mean(meanArray2)
# meanVar2 = np.std(meanArray2)

print "Mean 1: ",  mean1, " Variance: ", meanVar1
# print "Mean 2: ",  mean2, " Variance: ", meanVar2

