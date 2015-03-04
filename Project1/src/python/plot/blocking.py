import matplotlib.pyplot as plt
import numpy as np

def blocking(indata, outdata, blocksizes):
    infile = open(indata, 'r')
    outfile = open(outdata, 'w')

    array = np.loadtxt(infile)

    varray = []
    totarray = []
    smplarray = []

    totmean = np.mean(array)
    tmpmean_dif_tot=0
    tmpmean_dif_smpl=0
    
    for nb in blocksizes:
        tmplist = []
        
        for i in xrange(len(array)/nb):
            arr = array[nb*i:nb*(i+1)]
            tmp = np.mean(arr)
            tmplist.append(tmp)
            
            tmpmean_dif_tot += tmp-totmean
            for k in xrange(nb):
                tmpmean_dif_smpl += arr[k]-totmean
            
        tmpvar = np.var(tmplist)

        tmpmean_dif_tot = tmpmean_dif_tot/len(arr)
        tmpmean_dif_smpl = tmpmean_dif_smpl/len(array)

        varray.append(tmpvar)
        totarray.append(tmpmean_dif_tot)
        smplarray.append(tmpmean_dif_smpl)

        outfile.write(str(nb)+': '+' '+str(tmpmean_dif_tot)+' '\
                +str(tmpmean_dif_smpl)+' '+str(tmpvar)+'\n')

    infile.close()
    outfile.close()
    
    return blocksizes, varray, totarray, smplarray 

def blocks(magnitude):
    N = float(10**magnitude)
    blocksizes = []
    for i in range(magnitude+1):
        for j in range(magnitude+1):
            n = 5**i
            m = 2**j
            blocksizes.append(int(N/(m*n)))
    return blocksizes

x,y,tot,smpl = blocking('../../../res/energies.txt','blockvar.txt', np.sort(blocks(5)))

y = np.sort(y)
smpl = np.sort(smpl)
tot = np.sort(tot)


#plt.plot(x, np.sqrt(y))
plt.plot(x, tot)
plt.plot(x,smpl)
plt.legend(['tot','smpl'])
plt.show()
