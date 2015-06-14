import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/helium_04/"
resDir = "../../../../res/plot/helium_04/";

# Delete files that have the matching number as this program.
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if "02" in f]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass

beta = 0.3

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main_02.cpp"])
# Run the cpp code
call(["./a.out", str(beta)])
