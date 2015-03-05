from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

positions = open('../../../res/heliumWave1Wave2Density/wave1Positions.txt', 'r')
data = np.loadtxt(positions)

# Block size
nb = 50000 
N = len(data[0])/nb
particles = len(data)
r1Max = 3
r2Max = 3

bins = 100
r1 = np.linspace(0,r1Max,bins)
r2 = np.linspace(0,r2Max,bins)
y1 = np.zeros([bins,N])
y2 = np.zeros([bins,N])
y1sq = np.zeros([bins,N])
y2sq = np.zeros([bins,N])
std1 = np.zeros([bins])
std2 = np.zeros([bins])
y1Mean = np.zeros(bins)
y1MeanSq = np.zeros(bins)
y2Mean = np.zeros(bins)
y2MeanSq = np.zeros(bins)

# Create the N distribution plots
particle = 0
for i in xrange(N):
    for j in xrange(nb):
	value = data[particle][i*nb + j]
	if value > r1Max: continue
	index = int(value/r1Max*bins)
	y1[index][i] += 1;
    # y1 = y1/nb # Why don't I need this?

# # Create the squared sum array
# for i in xrange(bins):
#     for j in xrange(N):
# 	y1sq[i][j] += y1[i][j]*y1[i][j]

# Find the std for each bin
for i in xrange(bins):
    for j in xrange(N):
	y1Mean[i] += y1[i][j]
	y1MeanSq[i] += y1[i][j]*y1[i][j]
    y1Mean[i] = y1Mean[i]/N
    y1MeanSq[i] = y1MeanSq[i]/N
    std1[i] = np.sqrt(y1MeanSq[i] - y1Mean[i]*y1Mean[i])



y1Mean = y1Mean/nb
std1 = std1/nb

# Now the second particle
# Create the N distribution plots
particle = 1
for i in xrange(N):
    for j in xrange(nb):
	value = data[particle][i*nb + j]
	if value > r1Max: continue
	index = int(value/r1Max*bins)
	y2[index][i] += 1;

# Find the std for each bin
for i in xrange(bins):
    for j in xrange(N):
	y2Mean[i] += y2[i][j]
	y2MeanSq[i] += y2[i][j]*y2[i][j]
    y2Mean[i] = y2Mean[i]/N
    y2MeanSq[i] = y2MeanSq[i]/N
    std2[i] = np.sqrt(y2MeanSq[i] - y2Mean[i]*y2Mean[i])
y2Mean = y2Mean/nb
std2 = std2/nb


plt.plot(r1,y1Mean)
plt.plot(r1,y1Mean + std1)
plt.plot(r1,y1Mean - std1)

plt.show()

