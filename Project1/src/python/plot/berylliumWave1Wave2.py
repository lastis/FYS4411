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
        tmpmean = np.mean(tmplist)
        tmpvar = np.var(tmplist)
        
        print tmpvar

        # tmpmean_dif_tot = tmpmean_dif_tot
        # tmpmean_dif_smpl = tmpmean_dif_smpl

        varray.append(tmpvar)
        # totarray.append(tmpmean_dif_tot)
        # smplarray.append(tmpmean_dif_smpl)

        outfile.write(str(nb)+': '+str(tmpmean)+' '+str(tmpvar)+'\n')
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

x,y = blocking('../../../res/berylliumWave1Wave2/wave1Energies.txt', \
'../../../res/berylliumWave1Wave2/blockBerylliumWave1Energies.txt', np.sort(blocks(6)))

y = np.sort(y)

font = {'family' : 'serif',
        'size'   : 10}

plt.rc('font', **font)

fig, ax = plt.subplots()

ax.plot(x, np.sqrt(y), label='$E_0$ without Jastrow factor')

x,y = blocking('../../../res/berylliumWave1Wave2/wave2Energies.txt', \
'../../../res/berylliumWave1Wave2/blockBerylliumWave2Energies.txt', np.sort(blocks(6)))

y = np.sort(y)

ax.plot(x, np.sqrt(y), label='$E_0$ with Jastrow factor')

ax.set_title(r'Blocking Analysis of Ground State Energies as Function of STD with' '\n' 
    r'Generic Energy Calculation for Beryllium')

ax.legend(loc='best')
ax.set_xlabel(r'Blocksize', fontsize=14)
ax.set_ylabel(r'$\sigma$', fontsize=14)

ax.grid('on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#plt.savefig('../../../res/berylliumWave1Wave2/wave1EnergiesBlocking.eps')

plt.savefig('../../../res/berylliumWave1Wave2/wave12EnergiesBlocking.eps')
plt.show()
