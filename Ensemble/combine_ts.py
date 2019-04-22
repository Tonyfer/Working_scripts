lables = ['Left Caudate', 'Right Caudate', 'Left Putamen', 'Right Putamen', 'Left Thalamus', 'Right Thalamus', 'Left Accumbens', 'Right Accumbens', 'Left Pallidum', 'Right Pallidum']
import pickle as pkl
f = open('/data3/cdb/ytong/Ensemble/net_dic.pkl', 'rb')
net_dic = pkl.load(f)
f.close()


import os
import numpy as np
for folder in ['Ref/', 'Rep_50/']:
    sub_folder = '/data3/cdb/ytong/Ensemble/Sub/'+folder
    basc_folder = '/data3/cdb/ytong/Ensemble/BASC_444/'+folder
    list_sub = os.listdir(sub_folder)
    list_sub.sort()
    list_basc = os.listdir(basc_folder)
    list_basc.sort()

    for i in range(30):
        f = open(basc_folder+list_basc[i], 'rb')
        basc_ts = pkl.load(f)
        f.close()
        f = open(sub_folder+list_sub[i], 'rb')
        sub_ts = pkl.load(f)
        f.close()
        
        for j in range(1,8):
            part_basc = basc_ts[:, [i-1 for i in net_dic['network_'+str(j)]]]

            net = np.hstack((part_basc,sub_ts))
            
            f = open('/data3/cdb/ytong/Ensemble/yeo_network'+str(j)+'/'+folder+list_basc[i], 'wb')
            pkl.dump(net, f)
            f.close()
            







