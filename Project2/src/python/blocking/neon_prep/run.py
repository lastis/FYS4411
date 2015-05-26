import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/blocking/neon_prep/"
resDir = "../../../../res/blocking/neon_prep/";

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
alpha = 11.00

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main.cpp"])
# Run the cpp code
call(["./a.out", str(alpha), str(beta)])
