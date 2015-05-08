import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/berylliumAlphaQuantumNoCor/"
resDir = "../../../../res/plot/berylliumAlphaQuantumNoCor/";

fileName = "energies_mean.txt"
alphaFile = "alpha_array.txt"

# Delete the exsiting output file, if it exsists.
os.chdir(curDir)
os.chdir(resDir)
try :
    os.remove(fileName)
except OSError:
    pass

nCycles = 1e4
blockSize = nCycles/10
alphaArray = np.linspace(3.5,4.5,11)
# Save the alpha values.
np.savetxt(alphaFile,alphaArray)
# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main.cpp"])

os.chdir(curDir)
os.chdir(cppDir)
# Run the cpp code
for alpha in alphaArray:
    call(["./a.out", fileName, str(alpha),"1.0", str(nCycles), str(blockSize)])
