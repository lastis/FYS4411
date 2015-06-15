import os
from subprocess import call
from subprocess import Popen, PIPE
import numpy as np

curDir = os.getcwd();
# Adress to target directories from this directory. 
cppDir = "../../../cpp/plot/helium_06/"
resDir = "../../../../res/plot/helium_06/";


# Delete the exsiting .txt files
os.chdir(curDir)
os.chdir(resDir)
fileList = [f for f in os.listdir(".") if "01" in f]
try :
    for f in fileList:
        os.remove(f)
except OSError:
    pass

lengthA = 20
alphaArray = np.linspace(1,2.5,lengthA)
alphaFile = "alpha_array_01.txt"
np.savetxt(alphaFile,alphaArray)

lengthB = 20
betaArray = np.linspace(0.05,1.5,lengthB)
betaFile = "beta_array_01.txt"
np.savetxt(betaFile,betaArray)

energyArray = np.zeros([lengthA,lengthB])
energyFile = "energy_array_01.txt"

# Compile the code
os.chdir(curDir)
os.chdir(cppDir)
call(["make", "cpp_file = main_01.cpp"])

# Run the cpp code for helium with hydrogenlike wavefunctions
os.chdir(curDir)
os.chdir(cppDir)
counter = 0
for i, alpha in enumerate(alphaArray):
    for j, beta in enumerate(betaArray):
        # Run the cpp code
        p = Popen(["./a.out", str(alpha), str(beta)], stdin=PIPE,stdout=PIPE, stderr=PIPE)
        output, err = p.communicate(b"input data is passed to subprocess' stin")
        rc = p.returncode
        energy = float(output.split()[-1])
        # print "Beta:", beta, " Alpha: ", alpha, " Energy : ", energy
        counter = counter + 1
        print "cycle completed"
        # if counter == 10 :
        #     counter = 0
        #     print "10 Cycles completed."
        energyArray[i,j] = energy

os.chdir(curDir)
os.chdir(resDir)
np.savetxt(energyFile,energyArray)
