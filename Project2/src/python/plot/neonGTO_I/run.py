import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/neonGTO_I/"
resDir = "../../../../res/plot/neonGTO_I/";

betaFile = "beta_array.txt"

# Delete the exsiting .txt files
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if f.endswith(".txt")]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass

betaArray = np.linspace(0.1,4,11)

# Save the beta values.
np.savetxt(betaFile,betaArray)
# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main.cpp"])

os.chdir(curDir)
os.chdir(cppDir)
# Run the cpp code
for beta in betaArray:
    call(["./a.out", str(beta)])
