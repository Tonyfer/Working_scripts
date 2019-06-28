from nilearn import datasets
from nilearn import image
import nibabel as nib
import nympy as np
import pickle as pkl


basc = datasets.fetch_atlas_basc_multiscale_2015()
basc444 = nib.load(basc.scale444)

yeo = datasets.fetch_atlas_yeo_2011()
yeo7 = nib.load(yeo.thick_7)

template = datasets.load_mni152_template()
basc_re = image.resample_to_img(basc444, template, interpolation='nearest')
yeo_re = image.resample_to_img(yeo7, template, interpolation='nearest')


yeo_data = yeo_re.get_fdata()[:,:,:,0]
basc_data = basc_re.get_fdata()




inde = []
for i in range(1,445):
    li = list(np.unique(yeo_data[basc_data == i]))
    nu = [len(np.argwhere(yeo_data[basc_data == i] == j)) for j in li]
    inde.append(li[nu.index(max(nu))])

inde = np.array(inde)

net_dic = {}
for i in range(1,8):
    net_dic['network_'+str(i)] = list((np.argwhere(inde == i) + 1).reshape(-1))
    
    
    
f = open('net_dic.pkl', 'wb')
pkl.dump(net_dic,f)
f.close()









