# run as python re_label2.py ref target mask path_out smooth('Yes', 'No')



import numpy as np
from sklearn.preprocessing import OneHotEncoder
import nibabel as nib

def re_label2(ref, target):
    #on_hot encoding
    ref = ref.reshape(-1,1)
    target = target.reshape(-1,1)
    onehot_encoder = OneHotEncoder(sparse=False)
    ref_encoded = onehot_encoder.fit_transform(ref)
    target_encoded = onehot_encoder.fit_transform(target)
    dices = np.array([[dice(target_encoded[:,j], ref_encoded[:,i]) for i in range(ref_encoded.shape[1])] for j in range(target_encoded.shape[1])])
    
    pairs = []
    for i in range(dices.shape[0]): 
        pair = np.unravel_index(np.argmax(dices, axis=None), dices.shape)
        dices[pair[0],:] = -1
        dices[:,pair[1]] = -1
        pairs.append(pair)
    #loc = []
    #dic = {}
    #for j in range(target_encoded.shape[1]):
    #    ori_label = int(np.unique(target[target_encoded[:,j].astype(np.bool)]))
    #    dices = [dice(target_encoded[:,j], ref_encoded[:,i]) for i in range(ref_encoded.shape[1])]
    #    label = np.argmax(np.array(dices))
    #    label = int(np.unique(ref[ref_encoded[:,label].astype(np.bool)]))
    #    dic[str(ori_label)] = label
    #    target[target_encoded[:,j].astype(np.bool)] = label
    for i in pairs:
        target[target_encoded[:,i[0]].astype(np.bool)]= i[1]
    
    
    return target, pairs


def ndarray_to_vol(data_array, roi_mask_file, sample_file, filename):
    """
    Converts a numpy array to a nifti file given an roi mask
    Parameters
    ----------
    data_array : array_like
        A data array with the same column length and index alignment as the
        given roi_mask_file.  If data_array is two dimensional, first dimension
        is considered temporal dimension
    roi_mask_file : string
        Path of the roi_mask_file
    sample_file : string or list of strings
        Path of sample nifti file(s) to use for header of the output.
        If list, the first file is chosen.
    filename : string
        Name of output file
    Returns
    -------
    img_file : string
        Path of the nifti file output
    """

    import os
    import numpy as np
    import nibabel as nb

    roi_mask_file = nb.load(roi_mask_file).get_data().astype('bool')

    if data_array.ndim == 1:
        out_vol = np.zeros_like(roi_mask_file, dtype=data_array.dtype)
        out_vol[roi_mask_file] = data_array

    elif data_array.ndim == 2:
        list_roi_shape = list(roi_mask_file.shape[0:3])

        out_vol = np.zeros(
            list_roi_shape + [data_array.shape[1]],
            dtype=data_array.dtype
        )
        out_vol[roi_mask_file] = data_array

    else:
        raise ValueError(
            'data_array is %i dimensional, '
            'must be either 1 or 2 dimensional' % len(data_array.shape)
        )

    # TODO @AKI why not use header from ROI file?
    #           it should has the same affine
    if type(sample_file) is list:
        sample_file = sample_file[0]

    nii = nb.load(sample_file)

    img = nb.Nifti1Image(
        out_vol,
        header=nii.get_header(),
        affine=nii.get_affine()
    )

    img_file = os.path.join(os.getcwd(), filename)
    img.to_filename(img_file)

    return img_file, img



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
rel = re_label2(ref,target)[0]
rel = rel.reshape(shape)

mask = nib.load(sys.argv[3])

out = nib.Nifti1Image(rel,header=mask.get_header(),affine=mask.get_affine())
out.to_filename(sys.argv[4])


if sys.argv[5] == 'Yes':
    smooth(sys.argv[4], sys.argvs[4])
