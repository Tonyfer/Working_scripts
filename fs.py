# run with python fs.py [.txt file of fs] [data_type]




import pickle as pkl
from sklearn import preprocessing
import numpy as np

def normal_upper(array):
    upper = array[np.triu_indices(array.shape[0], k = 1)]
    return preprocessing.scale(upper)

def all_upper(sublist):
    return [normal_upper(i) for i in sublist]

import os
import sys
fs={}
with open(sys.argv[1],'r') as f:
    for line in f:
        fs[line.strip('\n').split(':')[0]] = line.strip('\n').split(':')[1].split(',')

        
        
        
cwd = '/data3/cdb/ytong/Ensemble/distance/'
for key, component in fs.items():
    value = []
    for i in component:
        f = open(cwd+sys.argv[2]+'_'+i+'.pkl','rb')  
        dis=pkl.load(f)  
        f.close()
        value.append(dis)
    output = np.sum(np.array([all_upper(i) for i in value]), axis=0)
    f1 = open('/data3/cdb/ytong/Ensemble/fs/'+sys.argv[2]+'_'+key+'.pkl', "wb")
    pkl.dump(output, f1)
    f1.close()
