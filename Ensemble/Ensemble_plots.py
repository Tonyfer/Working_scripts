import pickle as pkl
from sklearn import preprocessing
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import ptitprince as pt
import pandas as pd
from nilearn import plotting
from scipy import stats

def ICC(ratings):
    ns, nr = ratings.shape
    SStotal = ratings.var(ddof=1) * (ns * nr - 1)
    MSr = ratings.mean(axis=1).var(ddof=1) * nr
    MSw = (ratings.var(axis=1, ddof=1) / ns).sum()
    MSc = ratings.mean(axis=0).var(ddof=1) * ns
    MSe = (SStotal - MSr * (ns - 1) - MSc * (nr - 1)) / ((ns - 1) * (nr - 1))
    return (MSr - MSe)/(MSr + (nr - 1) * MSe)


def ICC_q(ratings):
    ns, nr = ratings.shape
    SStotal = ratings.var(ddof=1) * (ns * nr - 1)
    MSr = ratings.mean(axis=1).var(ddof=1) * nr
    MSw = (ratings.var(axis=1, ddof=1) / ns).sum()
    MSc = ratings.mean(axis=0).var(ddof=1) * ns
    MSe = (SStotal - MSr * (ns - 1) - MSc * (nr - 1)) / ((ns - 1) * (nr - 1))
    return MSr, MSe


length_list = [69, 59, 52, 46, 35, 57, 92]


for net in range(1,8):

    cwd = '/data3/cdb/ytong/Ensemble/yeo_network'+str(net)+'/fs/'
    plotd = '/data3/cdb/ytong/Ensemble/yeo_network'+str(net)+'/plot/'
#### Create matrix plot ####


    icc = []
    icc_msr = []
    icc_mse = []
    for k in ['Rep_50']:
        for i in ['fc','fs_5dis','fs_6dis']:
            ref = open(cwd+'Ref'+'_'+i+'.pkl','rb')
            fs_ref=pkl.load(ref)  
            ref.close()
            rep = open(cwd+k+'_'+i+'.pkl','rb')
            fs_rep=pkl.load(rep)  
            rep.close()
            icc += [ICC(np.transpose(np.array([fs_ref[:, i], fs_rep[:,i]]))) for i in range(fs_ref.shape[1])]
            icc_msr += [ICC_q(np.transpose(np.array([fs_ref[:, i], fs_rep[:,i]])))[0] for i in range(fs_ref.shape[1])]
            icc_mse += [ICC_q(np.transpose(np.array([fs_ref[:, i], fs_rep[:,i]])))[1] for i in range(fs_ref.shape[1])]
        

        
       
    fs = []
    for k in range(1):
        for i in ['fc','fs_5dis','fs_6dis']:
            fs+=[i for j in range(fs_ref.shape[1])]
        
    ref = []
    for k in ['Rep_50']:
        ref += [k for i in range(3*fs_ref.shape[1])]
    
    data = pd.DataFrame({'icc':icc, 'icc_msr':icc_msr, 'icc_mse':icc_mse, 'fs':fs, 'ref': ref})
        
    plt.figure(figsize=(20, 10))
    sns.pointplot(x="ref", y="icc_msr", data=data, hue= 'fs', dodge=0.53, join=False, palette="dark",markers="d", scale=.75, ci='sd',capsize = 0.07)
    sns.stripplot(x="ref", y="icc_msr", data=data, hue= 'fs', size = 3, dodge=0.45, alpha = 0.05).set_title('Edge-wise ICC MSr')
    pt.half_violinplot(x="ref", y="icc_msr", data=data, hue= 'fs',scale = "area",inner = None, offset = 0.03, saturation=0.5)
    plt.legend(ncol=2)       
    plt.savefig(plotd+'icc_msr.png')######
    plt.close()

    plt.figure(figsize=(20, 10))
    sns.pointplot(x="ref", y="icc_mse", data=data, hue= 'fs', dodge=0.53, join=False, palette="dark",markers="d", scale=.75, ci='sd',capsize = 0.07)
    sns.stripplot(x="ref", y="icc_mse", data=data, hue= 'fs', size = 3, dodge=0.45, alpha = 0.05).set_title('Edge-wise ICC MSe')
    pt.half_violinplot(x="ref", y="icc_mse", data=data, hue= 'fs',scale = "area",inner = None, offset = 0.03, saturation=0.5)

    t1,p1 = stats.ttest_ind(icc_mse[0:int(len(icc_mse)/3)],icc_mse[int(len(icc_mse)/3):int(len(icc_mse)/3*2)], nan_policy ='omit', equal_var=False)
    t2,p2 = stats.ttest_ind(icc_mse[0:int(len(icc_mse)/3)],icc_mse[int(len(icc_mse)/3*2):], nan_policy ='omit', equal_var=False)

    plt.text(-0.15,0,'T: '+str(round(t1,5))+'\n'+'P: '+str(round(p1,5)),fontsize=18)
    plt.text(0.15,0,'T: '+str(round(t2,5))+'\n'+'P: '+str(round(p2,5)),fontsize=18)
    plt.legend(ncol=2)
    plt.savefig(plotd+'icc_mse.png')######
    plt.close()
        
      
        
        

        
    data_1 = icc[0:fs_ref.shape[1]]
    out_1 = np.zeros((length_list[net-1], length_list[net-1]))
    inds = np.triu_indices(len(out_1), k = 1)
    out_1[inds] = data_1

    plotting.plot_matrix(out_1, figure=(10, 8),reorder=False, vmax=1, vmin=-1)
    plt.title('ICC')
    plt.savefig(plotd+'icc_new.png')#####
    plt.close()

    data_2 = icc[(fs_ref.shape[1]):(fs_ref.shape[1]*2)]
    out_2 = np.zeros((length_list[net-1], length_list[net-1]))
    inds = np.triu_indices(len(out_1), k = 1)
    out_2[inds] = data_2

    plotting.plot_matrix(out_2-out_1, figure=(10, 8),reorder=False,vmax=0.3, vmin=-0.3)
    plt.title('ICC difference between fc & 5_dis')
    plt.savefig(plotd+'mat_dif_5.png')#####
    plt.close()

    data_3 = icc[(fs_ref.shape[1]*2):]
    out_3 = np.zeros((length_list[net-1], length_list[net-1]))
    inds = np.triu_indices(len(out_1), k = 1)
    out_3[inds] = data_3

    plotting.plot_matrix(out_3-out_1, figure=(10, 8),reorder=False,vmax=0.3, vmin=-0.3)
    plt.title('ICC difference between fc & 6_dis')
    plt.savefig(plotd+'mat_dif_6.png')#####
    plt.close()



    diff_1 = np.array(data_1)-np.array(data_2)
    diff_2 = np.array(data_1)-np.array(data_3)

    sns.distplot(diff_1, bins=20, kde=True).set(xlim=(-0.15, 0.15))
    plt.title('hist_of_5_dis')
    plt.savefig(plotd+'hist_dif_5.png')#####
    plt.close()

    sns.distplot(diff_2, bins=20, kde=True).set(xlim=(-0.15, 0.15))
    plt.title('hist_of_6_dis')
    plt.savefig(plotd+'hist_dif_6.png')#####
    plt.close()





















'''

#### Create Sub-wise Absolute Difference Plot ####

cwd = '/data3/cdb/ytong/Ensemble/BASC_020/fs/'
value = []
for k in ['Rep_50']:
    for i in ['fc','fs_5dis','fs_6dis']:
        ref = open(cwd+'Ref'+'_'+i+'.pkl','rb')
        fs_ref=pkl.load(ref)  
        ref.close()
        rep = open(cwd+k+'_'+i+'.pkl','rb')
        fs_rep=pkl.load(rep)  
        rep.close()
        value += list(np.mean(abs(fs_rep-fs_ref), axis = 1))
        
fs = []
for k in range(1):
    for i in ['fc','fs_5dis','fs_6dis']:
        fs+=[i for j in range(fs_ref.shape[0])]
        
ref = []
for k in ['Rep_50']:
    ref += [k for i in range(3*fs_ref.shape[0])]
    
data = pd.DataFrame({'value':value, 'fs':fs, 'ref': ref})

t1,p1 = stats.ttest_ind(value[0:int(len(value)/3)],value[int(len(value)/3):int(len(value)/3*2)], nan_policy ='omit', equal_var=False)

t2,p2 = stats.ttest_ind(value[0:int(len(value)/3)],value[int(len(value)/3*2):], nan_policy ='omit', equal_var=False)

plt.figure(figsize=(20, 10))
sns.pointplot(x="ref", y="value", data=data, hue= 'fs', dodge=0.53, join=False, palette="dark", markers="d", scale=.75, ci='sd',capsize = 0.07)
sns.stripplot(x="ref", y="value", data=data, hue= 'fs', size = 3, dodge=0.45, alpha = 0.5).set_title('Sub-wise Absolute Difference 020')
#sns.violinplot(x="fs", y="value", data=data, hue= 'ref', inner =  None, split =True) 
pt.half_violinplot(x="ref", y="value", data=data, hue= 'fs', scale = "area",inner = None, offset = 0.03, saturation=0.5)

plt.text(-0.15,0.15,'T: '+str(round(t1,5))+'\n'+'P: '+str(round(p1,5)))
plt.text(0.15,0.15,'T: '+str(round(t2,5))+'\n'+'P: '+str(round(p2,5)))

plt.legend(ncol=2)





#plt.savefig('test.png')


'''
