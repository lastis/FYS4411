import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/berylliumAlpha/"
resDir = "../../../../res/plot/berylliumAlpha/";

fileName = "energies_mean.txt"
alphaFile = "alpha_array.txt"

# Delete the exsiting output file. 
os.chdir(curDir)
os.chdir(resDir)
try :
    os.remove(fileName)
except OSError:
    pass
# Save the alpha values.
alphaArray = np.linspace(2,5,5)
np.savetxt(alphaFile,alphaArray)
# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main.cpp"])

os.chdir(curDir)
os.chdir(cppDir)
# Run the cpp code
for alpha in alphaArray:
    call(["./a.out", fileName, str(alpha),"0.8", "1000"])
