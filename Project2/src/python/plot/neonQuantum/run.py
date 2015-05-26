import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/neonQuantum/"
resDir = "../../../../res/plot/neonQuantum/";

alphaFile = "alpha_array.txt"

# Delete the exsiting .txt files
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if f.endswith(".txt")]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass

beta = 0.21
alphaArray = np.linspace(9.0,11.0,11)


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
    call(["./a.out", str(alpha), str(beta)])
