import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/heliumBeta/"
resDir = "../../../../res/plot/heliumBeta/";

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

alpha = 1.66
betaArray = np.linspace(0.1,2.0,21)


# Save the alpha values.
np.savetxt(betaFile,betaArray)
# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make"])

os.chdir(curDir)
os.chdir(cppDir)
# Run the cpp code
for beta in betaArray:
    call(["./a.out", str(alpha), str(beta)])
