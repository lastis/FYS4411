import numpy as np

x = np.random.uniform(0,1,10e5)

outfile = open('test.txt','w')

for i in x:
    outfile.write(str(i)+'\n')

outfile.close()
