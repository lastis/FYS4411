import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/beryllium_01/"
resDir = "../../../../res/plot/beryllium_01/";


# Delete the exsiting .txt files
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if "02" in f]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass

betaArray = np.linspace(0.001,0.301,21)

# Save the beta values.
betaFile = "beta_array_02.txt"
np.savetxt(betaFile,betaArray)

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make", "cpp_file = main_02.cpp"])

# Run the cpp code for helium with gto wave functions. 
os.chdir(curDir)
os.chdir(cppDir)
for beta in betaArray:
    call(["./a.out", str(beta)])
