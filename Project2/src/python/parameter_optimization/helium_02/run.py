import os
from subprocess import call
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/parameter_optimization/helium_02/"
resDir = "../../../../res/parameter_optimization/helium_02/";


# Find optimal values for helium with STO's, not using importance sampling. 
# Use CG method, with static N, but n_b (block size) is chosen from the
# blocking plot to find an optimal value to find E_0. 

# Variables for finding the numerical approximation to
# the energy gradient from the paramters alpha and beta. 
def f(x,*args):
    beta = x[0]
    if (beta < 0):
        beta = -beta
    # Run the cpp code
    p = Popen(["./a.out",str(beta)], stdin=PIPE,stdout=PIPE, stderr=PIPE)
    output, err = p.communicate(b"input data is passed to subprocess' stin")
    rc = p.returncode
    print "Beta:", beta
    print output
    # Get the third last word of the output which is the energy.
    energy = float(output.split()[-2])
    return energy

def gradf(x,*args):
    beta = x[0]
    if (beta < 0):
        beta = -beta
    # Run the cpp code
    p = Popen(["./a.out", str(beta)], stdin=PIPE,stdout=PIPE, stderr=PIPE)
    output, err = p.communicate(b"input data is passed to subprocess' stin")
    rc = p.returncode
    print "Beta:", beta
    print output
    # Get the two last word of the output which is the the gradient
    output = output.split()
    dE_dBeta = float(output[-1])
    return np.asarray(dE_dBeta)

# START OF PROGRAM

# Delete the exsiting .txt files
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if f.endswith(".txt")]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make"])

# Set the initial values for the minimization function.
x0 = np.asarray(0.35)

# Run the minimization code. 
res = optimize.fmin_cg(f,x0,fprime=gradf,disp=True)
print res

# Output results to a file.
os.chdir(resDir)
myfile = open("output.txt", "w")
myfile.write(res)
myfile.close()

