import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/helium_04/"
resDir = "../../../../res/plot/helium_04/";


# Delete the exsiting .txt files
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if "01" in f]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass


alpha = 1.79
beta = 0.4125

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make", "cpp_file = main_01.cpp"])

# Run the cpp code for helium with hydrogenlike wavefunctions
os.chdir(curDir)
os.chdir(cppDir)
call(["./a.out", str(alpha), str(beta)])
