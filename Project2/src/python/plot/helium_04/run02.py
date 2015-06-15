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
fileList = [f for f in os.listdir(".") if "02" in f]
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
