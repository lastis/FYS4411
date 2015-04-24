from subprocess import call
import os

curDir = os.getcwd();
resDir = "../../../../res/";
cppDir = "../../../cpp/plot/"

os.chdir(cppDir)
call(["make","cpp_file=berylliumAlphas.cpp"])
call(["./a.out", "4","0.8", "50"])

