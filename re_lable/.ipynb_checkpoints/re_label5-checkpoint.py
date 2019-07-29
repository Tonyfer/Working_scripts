# run as python re_label2.py ref target mask path_out smooth('Yes', 'No')



import numpy as np
from sklearn.preprocessing import OneHotEncoder
import nibabel as nib
import sys
import scipy

def dice(im1, im2, empty_score=1.0):
   
    im1 = np.asarray(im1).astype(np.bool)
    im2 = np.asarray(im2).astype(np.bool)

    if im1.shape != im2.shape:
        raise ValueError("Shape mismatch: im1 and im2 must have the same shape.")

    im_sum = im1.sum() + im2.sum()
    if im_sum == 0:
        return empty_score

    # Compute Dice coefficient
    intersection = np.logical_and(im1, im2)

    return 2*intersection.sum() / im_sum



def re_label4(ref, target):
    from itertools import permutations
    #on_hot encoding
    ref = ref.reshape(-1,1)
    target = target.reshape(-1,1)
    onehot_encoder = OneHotEncoder(sparse=False)
    ref_encoded = onehot_encoder.fit_transform(ref)
    target_encoded = onehot_encoder.fit_transform(target)
    dices = np.array([[dice(target_encoded[:,j], ref_encoded[:,i]) for i in range(ref_encoded.shape[1])] for j in range(target_encoded.shape[1])])
    count = [np.sum(target == i) for i in range(dices.shape[0])]
    ind = np.argsort(count)[::-1]
    
    
    pairs = []
    for i in ind:
        match = np.argmax(dices[i, :])
        dices[i,:] = -1
        dices[:, match] = -1
        pairs.append((i,match))
    
    
    for i in pairs:
        target[target_encoded[:,i[0]].astype(np.bool)]= i[1]
    
    
    return target, pairs








def smooth(path_in, path_out):
    ref = nib.load(path_in)
    header = ref.get_header()
    affine = ref.get_affine()
    ref = ref.get_fdata()
    ooo = np.zeros(ref.shape)
    for i in range(1, ref.shape[0]-1):
        for j in range(1, ref.shape[1]-1):
            for k in range(1, ref.shape[2]-1):
                if ref[i,j,k] == 0:
                    ooo[i,j,k] = 0
                else:
                    cube = ref[(i-1):(i+2), (j-1):(j+2), (k-1):(k+2)]
                    if np.sum(cube!=0) == 0:
                        q = 0
                    else:
                
            #un = np.array([np.sum(cube == d) for d in list(np.unique(cube))])
                        q = int(scipy.stats.mode(cube[cube>0],axis=None)[0])
                    ooo[i,j,k] = q
    
    
    
    out = nib.Nifti1Image(ooo,header=header,affine=affine)
    out.to_filename(path_out)
    
    

ref = nib.load(sys.argv[1])
ref = ref.get_fdata()
shape = ref.shape
ref = ref.reshape(-1,1)
target = nib.load(sys.argv[2])
target = target.get_fdata().reshape(-1,1)
rel = re_label4(ref,target)[0]
rel = rel.reshape(shape)

mask = nib.load(sys.argv[3])

out = nib.Nifti1Image(rel,header=mask.get_header(),affine=mask.get_affine())
out.to_filename(sys.argv[4])


if sys.argv[5] == 'Yes':
    smooth(sys.argv[4], sys.argv[4])