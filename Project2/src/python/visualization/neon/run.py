import os
from subprocess import call
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/visualization/neon/"
resDir = "../../../../res/visualization/neon/";
fileName = "density3D.txt"

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main.cpp"])

alpha = 9.82
beta = 0.18

call(["./a.out", str(alpha), str(beta)])




