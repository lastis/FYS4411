import os
from subprocess import call
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/parameter_optimization/helium_02/"
resDir = "../../../../res/parameter_optimization/helium_02"


# Find optimal values for helium with STO's, not using importance sampling. 
# Use CG method, with static N, but n_b (block size) is chosen from the
# blocking plot to find an optimal value to find E_0. 

# Variables for finding the numerical approximation to
# the energy gradient from the paramters alpha and beta. 

def func(beta):
    # Run the cpp code
    p = Popen(["./a.out", str(beta)], stdin=PIPE,stdout=PIPE, stderr=PIPE)
    output, err = p.communicate(b"input data is passed to subprocess' stin")
    rc = p.returncode
    # print output
    # Get the third last word of the output which is the energy.
    energy = float(output.split()[-2])
    return energy

# START OF PROGRAM

# Delete files that have the matching number as this program.
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
call(["make"])

# Set the initial values for the minimization function.
betaOld = 0.252
dBeta = 0.1
dBetaThr = 0.001

# Run the minimization code. 
energyOld = func(betaOld)
canImprove = True
while canImprove:
    print "Energy : ", energyOld,  " Beta: ", betaOld
    if dBeta < dBetaThr : 
        break

    betaNew = betaOld - dBeta
    energyNew = func(betaNew)
    if energyNew < energyOld :
        betaOld = betaNew
        energyOld = energyNew
        print "beta -"
        continue;

    betaNew = betaOld + dBeta
    energyNew = func(betaNew)
    if energyNew < energyOld :
        betaOld = betaNew
        energyOld = energyNew
        print "beta + "
        continue;

    dBeta = dBeta / 10
    print "Cycle completed"


# # Output results to a file.
# os.chdir(resDir)
# myfile = open("output.txt", "w")
# myfile.write(res)
# myfile.close()

