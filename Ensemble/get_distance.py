#run with python get_distance.py data_directory distance_type time_lag



import sklearn.metrics
import numpy as np

#def distance(data_type, req):
 #   dis = [sklearn.metrics.pairwise_distances(np.transpose(ts[0:(ts.shape[0]-i), :]), #Y=np.transpose(ts[i:,:]), metric=req) for i in range(6)]
#    return dis

import sys
import os
import pickle as pkl

cwd = sys.argv[1]
files = os.listdir(cwd)



ts_d = []
for file in files:
    f=open(cwd+file,'rb')  
    ts=pkl.load(f)  
    f.close()
    dis = sklearn.metrics.pairwise_distances(np.transpose(ts[0:(ts.shape[0]-int(sys.argv[3])), :]), Y=np.transpose(ts[int(sys.argv[3]):,:]), metric=sys.argv[2])
    ts_d.append(dis)
    
    
    
f1 = open('/data3/cdb/ytong/Ensemble/distance/'+sys.argv[1].split('/')[-2]+'_'+sys.argv[2]+'_t' + sys.argv[3]+'.pkl', "wb")
pkl.dump(ts_d, f1)
f1.close()
    










