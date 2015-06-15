import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/neon_02/"
resDir = "../../../../res/plot/neon_02/";

# Delete files that have the matching number as this program.
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if "01" in f]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass

beta = 0.1
alpha = 10.0

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main_01.cpp"])
# Run the cpp code
call(["./a.out", str(alpha), str(beta)])
