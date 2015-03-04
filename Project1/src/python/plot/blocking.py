import matplotlib.pyplot as plt
import numpy as np

def blocking(indata, outdata, blocksizes):
    infile = open(indata, 'r')
    outfile = open(outdata, 'w')

    array = np.loadtxt(infile)
    varray = []
    
    # totarray = []
    # smplarray = []

    # totmean = np.mean(array)
    # tmpmean_dif_tot=0
    # tmpmean_dif_smpl=0
    
    for nb in blocksizes:
        tmplist = []
        N = len(array)
        lim = N/nb
        print lim
        for i in xrange(lim): 
            arr = array[nb*i:nb*(i+1)]
            tmp = np.mean(arr)
            tmplist.append(tmp)
            # tmpmean_dif_tot += tmp-totmean
            
            # for k in xrange(nb):
            #     tmpmean_dif_smpl += arr[k]-totmean    
        
        tmpvar = np.var(tmplist)
        
        print tmpvar

        # tmpmean_dif_tot = tmpmean_dif_tot
        # tmpmean_dif_smpl = tmpmean_dif_smpl

        varray.append(tmpvar)
        # totarray.append(tmpmean_dif_tot)
        # smplarray.append(tmpmean_dif_smpl)

        outfile.write(str(nb)+': '+' '+str(tmpvar)+'\n')
        # str(tmpmean_dif_tot)+' '+str(tmpmean_dif_smpl)+')

    infile.close()
    outfile.close()
    
    return blocksizes, varray # , totarray, smplarray 

def blocks(magnitude):
    N = float(10**magnitude)
    blocksizes = []
    for i in range(magnitude+1):
        for j in range(magnitude+1):
            n = 5**i
            m = 2**j
            blocksizes.append(int(N/(m*n)))
    return blocksizes

def blocks2(m):
    blocksizes = []
    i = 0
    while i < 10**m:
        i += 10000
        blocksizes.append(i)
    return blocksizes

print blocks(6)

x,y = blocking('../../../res/heliumWave1Wave2/wave1Energies.txt','blockVarwave1Energies.txt', np.sort(blocks(6)))

y = np.sort(y)

plt.plot(x, np.sqrt(y))
plt.title(r'Blocking Analysis of Ground State Energies as Function of STD with' '\n' 
    r'Generic Energy Calculation and Trial Wavefunction $\psi_{T1}$')
plt.xlabel(r'Blocksize', fontsize=14)
plt.ylabel(r'$\sigma$', fontsize=14)
plt.grid('on')

plt.savefig('../../../res/heliumWave1Wave2/wave1EnergiesBlocking.eps')
plt.show()
plt.figure()
