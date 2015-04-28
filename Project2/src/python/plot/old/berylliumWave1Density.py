from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

positions = open('../../../res/berylliumWave1Wave2Density/wave1Positions.txt', 'r')
data = np.loadtxt(positions)

# Block size
nb = 50000 
N = len(data[0])/nb
particles = len(data)
bins = 100

r1Max = 3
r2Max = 3
r3Max = 10
r4Max = 10

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

r3 = np.linspace(0,r3Max,bins)
r4 = np.linspace(0,r4Max,bins)
y3 = np.zeros([bins,N])
y4 = np.zeros([bins,N])
y3sq = np.zeros([bins,N])
y4sq = np.zeros([bins,N])
std3 = np.zeros([bins])
std4 = np.zeros([bins])
y3Mean = np.zeros(bins)
y3MeanSq = np.zeros(bins)
y4Mean = np.zeros(bins)
y4MeanSq = np.zeros(bins)


# Create the N distribution plots
particle = 0
for i in xrange(N):
    for j in xrange(nb):
	value = data[particle][i*nb + j]
	if value >= r1Max: continue
	index = int(value/r1Max*bins)
	y1[index][i] += 1;

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

print "particle 1 finished"

# Now the second particle
# Create the N distribution plots
particle = 1
for i in xrange(N):
    for j in xrange(nb):
	value = data[particle][i*nb + j]
	if value >= r2Max: continue
	index = int(value/r2Max*bins)
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

print "particle 2 finished"

# Now the third particle
# Create the N distribution plots
particle = 2
for i in xrange(N):
    for j in xrange(nb):
	value = data[particle][i*nb + j]
	if value >= r3Max: continue
	index = int(value/r3Max*bins)
	y3[index][i] += 1;

# Find the std for each bin
for i in xrange(bins):
    for j in xrange(N):
	y3Mean[i] += y3[i][j]
	y3MeanSq[i] += y3[i][j]*y3[i][j]
    y3Mean[i] = y3Mean[i]/N
    y3MeanSq[i] = y3MeanSq[i]/N
    std3[i] = np.sqrt(y3MeanSq[i] - y3Mean[i]*y3Mean[i])
y3Mean = y3Mean/nb
std3 = std3/nb

print "particle 3 finished"
# Now the fourth particle
# Create the N distribution plots
particle = 3
for i in xrange(N):
    for j in xrange(nb):
	value = data[particle][i*nb + j]
	if value >= r4Max: continue
	index = int(value/r4Max*bins)
	y4[index][i] += 1;

# Find the std for each bin
for i in xrange(bins):
    for j in xrange(N):
	y4Mean[i] += y4[i][j]
	y4MeanSq[i] += y4[i][j]*y4[i][j]
    y4Mean[i] = y4Mean[i]/N
    y4MeanSq[i] = y4MeanSq[i]/N
    std4[i] = np.sqrt(y4MeanSq[i] - y4Mean[i]*y4Mean[i])
y4Mean = y4Mean/nb
std4 = std4/nb
print "particle 4 finished"


plt.plot(r1,y1Mean)
plt.plot(r1,y2Mean)
plt.plot(r1,y3Mean)
plt.plot(r1,y4Mean)

plt.show()

