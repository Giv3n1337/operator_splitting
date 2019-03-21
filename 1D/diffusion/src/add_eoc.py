import os
import sys
import numpy as np

dataDir = str(sys.argv[1])
if dataDir == '.':
    dataDir = os.getcwd()
if dataDir[-1] != '/':
    dataDir = dataDir + '/'

fileNames = os.listdir(dataDir)

for fileName in fileNames:
    data = np.loadtxt(dataDir+fileName, dtype=float)
    eoc = np.zeros(data.shape[0], dtype=float)

    for i in range(1,data.shape[0]):
        eoc[i] = np.log(data[i,1] / data[i-1,1]) / np.log(data[i-1,0] / data[i,0])
            
    data = np.concatenate((data, eoc[:,np.newaxis]), axis=1)
    
    np.savetxt(dataDir+fileName, data, fmt=('%4i', '%.16f', '%1.2f'), delimiter = 2*'\t')
        
