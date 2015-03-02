import matplotlib.pyplot as plt
import numpy as np

def blocking(indata, outdata, blocksizes):
    infile = open(indata, 'r')
    outfile = open(outdata, 'w')

    array = np.loadtxt(infile)

    meanvector = []

    for nb in blocksizes:
        tmplist = []
        
        for i in range(len(array)/nb):
            tmp = np.mean(array[nb*i:nb*(i+1)])
            tmplist.append(tmp)
            #outfile.write(str(tmp)+' ')

        tmpvar = np.var(tmplist)
        tmpmean = np.mean(tmplist)

        outfile.write(str(nb)+': '+' '+str(tmpmean)+' '+str(tmpvar)+'\n')
        meanvector.append(tmplist)

    infile.close()
    outfile.close()
    
    return 0

def blocks(magnitude):
    N = float(10**magnitude)
    blocksizes = []
    for i in range(magnitude+1):
        for j in range(magnitude+1):
            n = 5**i
            m = 2**j
            blocksizes.append(int(N/(m*n)))
    return blocksizes

blocking('test.txt','testdump.txt', np.sort(blocks(5)))
