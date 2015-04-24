from subprocess import call
import os

curDir = os.getcwd();
resDir = "../../../../res/";
cppDir = "../../../cpp/plot/"

os.chdir(cppDir)
call(["make","cpp_file=berylliumAlpha.cpp"])
call(["./a.out"])

