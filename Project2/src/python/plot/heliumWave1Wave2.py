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

fig1, ax1 = plt.subplots()

## Closed-form
x,y = blocking('../../../res/heliumWave2ClosedImportanceGeneric/closedEnergies.txt','blockVarwave1Energies.txt', np.sort(blocks(6)))
y = np.sort(y)
ax1.plot(x, np.sqrt(y),label='Closed-form energy')


## Closed-form with Imp
x,y = blocking('../../../res/heliumWave2ClosedImportanceGeneric/closedImportanceEnergies.txt',
	'blockVarwave2ciEnergies.txt', np.sort(blocks(6)))
y = np.sort(y)
ax1.plot(x, np.sqrt(y),label='Closed-form energy w/imp. sampling')


## Generic energy
x,y = blocking('../../../res/heliumWave2ClosedImportanceGeneric/genericEnergies.txt',
	'blockVarwave2genEnergies.txt', np.sort(blocks(6)))
y = np.sort(y)
ax1.plot(x, np.sqrt(y),label='Generic energy')


ax1.set_title(r'Blocking Analysis of Ground State Energies as Function of STD with' '\n' 
    r'Trial Wavefunction $\psi_{T2}$')
ax1.set_xlabel(r'Blocksize', fontsize=14)
ax1.set_ylabel(r'$\sigma$', fontsize=14)
ax1.legend(loc='best')

ax1.grid('on')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

plt.savefig('../../../res/heliumWave2ClosedImportanceGeneric/heliumWave2BlockingClosedImpGen.pdf')
plt.show()
plt.figure()


