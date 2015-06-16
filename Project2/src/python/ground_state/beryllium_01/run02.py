import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/ground_state/beryllium_01/"
resDir = "../../../../res/ground_state/beryllium_01/";


# Delete the exsiting .txt files
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if "energies_02" in f]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass


# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make", "cpp_file = main_02.cpp"])

# Run the cpp code for helium with hydrogenlike wavefunctions
os.chdir(curDir)
os.chdir(cppDir)
call(["./a.out"])
