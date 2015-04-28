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



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.plot(alpha, energy)

ax1.set_title(r'Ground State Energies as Function of $\alpha$ for Helium with' '\n' 
    r'Generic Energy Calculation and Trial Wavefunction $\psi_{T1}$')
ax1.set_xlabel(r'$\alpha$', fontsize=14)
ax1.set_ylabel(r'$E_0$', fontsize=14)

ax1.grid('on')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

plt.savefig('../../../res/heliumWave1Alpha/heliumWave1Alpha_alpha.eps')
plt.show()
plt.figure()



fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

ax2.plot(alpha,variance)

ax2.set_title(r'Variance as Function of $\alpha$ for Helium with' '\n' 
    r'Generic Energy Calculation and Trial Wavefunction $\psi_{T1}$')
ax2.set_xlabel(r'$\alpha$', fontsize=14)
ax2.set_ylabel(r'$Variance$', fontsize=14)

ax2.grid('on')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.savefig('../../../res/heliumWave1Alpha/heliumWave1Alpha_variance.eps')
plt.show()
plt.figure()
