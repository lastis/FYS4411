import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/berylliumBlocking/"
resDir = "../../../../res/plot/berylliumBlocking/";

fileName = "std_vs_blocksize.txt"

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main.cpp"])
# Run the cpp code
alpha = 3.5
call(["./a.out", fileName, str(alpha),"0.8"])
