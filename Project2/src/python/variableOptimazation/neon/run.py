import os
from subprocess import call
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/variableOptimazation/neon/"
resDir = "../../../../res/variableOptimazation/neon/";
fileName = "energies_mean.txt"
alphaFile = "alpha_array.txt"

ah = 0.01
ahInv = 1/ah
bh = 0.001
bhInv = 1/bh

def f(x,*args):
    alpha, beta = x
    if beta < 0 : beta = -beta
    # Run the cpp code
    os.chdir(curDir)
    os.chdir(cppDir)
    call(["./a.out", str(alpha), str(beta)])
    # p = Popen(["./a.out", str(alpha), str(beta)], stdin=PIPE,stdout=PIPE, stderr=PIPE)
    # output, err = p.communicate(b"input data is passed to subprocess' stin")
    # rc = p.returncode
    # energy = float(output)

    # Read output from file
    os.chdir(curDir)
    os.chdir(resDir)
    energyFile    = open(fileName, 'r')
    energyArray   = np.loadtxt(energyFile)
    energyFile.close()
    energy = np.mean(energyArray)
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
alpha = 9.82
beta = 0.18

x0 = np.asarray((alpha,beta))
args = 0
res = optimize.fmin_ncg(f,x0,fprime=gradf,disp=True)
print res



