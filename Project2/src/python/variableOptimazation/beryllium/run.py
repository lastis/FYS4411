import os
from subprocess import call
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/variableOptimazation/beryllium/"
resDir = "../../../../res/variableOptimazation/beryllium/";
fileName = "energies_mean.txt"
alphaFile = "alpha_array.txt"

ah = 0.1
ahInv = 10
bh = 0.01
bhInv = 100

def f(x,*args):
    alpha, beta = x
    # Run the cpp code
    p = Popen(["./a.out", str(alpha), str(beta)], stdin=PIPE,stdout=PIPE, stderr=PIPE)
    output, err = p.communicate(b"input data is passed to subprocess' stin")
    rc = p.returncode
    energy = float(output)
    print energy, alpha, beta
    return energy

def gradf(x,*args):
    alpha, beta = x
    galpha = (f((alpha+ah,beta),0) - f((alpha-ah,beta),0)) * 0.5*ahInv
    gbeta = (f((alpha,beta+bh),0) - f((alpha,beta-bh),0)) * 0.5*bhInv
    return np.asarray((galpha,gbeta))

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make","cpp_file=main.cpp"])

x0 = np.asarray((3.6,1.57))
args = 0
# res = optimize.fmin_cg(f,x0)
res = optimize.fmin_cg(f,x0,fprime=gradf,disp=True)
print res

os.chdir(resDir)
myfile = open("output.txt", "w")
myfile.write(res)
myfile.close()

