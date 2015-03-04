import matplotlib.pyplot as plt
import numpy as np

infile = open('../../../res/heliumWave1Alpha/alpha_energy_variance.txt', 'r')

heliumWave1Alpha = np.loadtxt(infile)

alpha = heliumWave1Alpha[0]

energy = heliumWave1Alpha[1]

variance = heliumWave1Alpha[2]

infile.close()

font = {'family' : 'serif',
        'size'   : 10}

plt.rc('font', **font)

plt.plot(alpha, energy)
plt.title(r'Ground State Energies as Function of $\alpha$ for Helium with' '\n' 
    r'Generic Energy Calculation and Trial Wavefunction $\psi_{T1}$')
plt.xlabel(r'$\alpha$', fontsize=14)
plt.ylabel(r'$E_0$', fontsize=14)
plt.grid('on')

plt.savefig('../../../res/heliumWave1Alpha/heliumWave1Alpha_alpha.eps')
plt.show()
plt.figure()

plt.plot(alpha, variance)
plt.title(r'Variance as Function of $\alpha$ for Helium with' '\n' 
    r'Generic Energy Calculation and Trial Wavefunction $\psi_{T1}$')
plt.xlabel(r'$\alpha$', fontsize=14)
plt.ylabel(r'$E_0$', fontsize=14)
plt.grid('on')

plt.savefig('../../../res/heliumWave1Alpha/heliumWave1Alpha_variance.eps')
plt.show()
plt.figure()
