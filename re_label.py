import numpy as np
from sklearn.preprocessing import OneHotEncoder


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




def re_label(ref, target):
    #on_hot encoding
    ref = ref.reshape(-1,1)
    target = target.reshape(-1,1)
    onehot_encoder = OneHotEncoder(sparse=False)
    ref_encoded = onehot_encoder.fit_transform(ref)
    target_encoded = onehot_encoder.fit_transform(target)
    loc = []
    dic = {}
    
    #compare the dice coefficients of different labels
    for j in range(target_encoded.shape[1]):
        ori_label = int(np.unique(target[target_encoded[:,j].astype(np.bool)]))
        dices = [dice(target_encoded[:,j], ref_encoded[:,i]) for i in range(ref_encoded.shape[1])]
        label = np.argmax(np.array(dices))
        label = int(np.unique(ref[ref_encoded[:,label].astype(np.bool)]))
        dic[str(ori_label)] = label
        target[target_encoded[:,j].astype(np.bool)] = label
        
    # return re_labeled array and the dictionary of label mapping 
    return target, dic


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






"""
usage:

ref = np.lpoad(ref_path)
target = np.load(target_path)

relabeled, label_dic = re_label(ref, target)
file,_= ndarray_to_vol(relabeled, 'BG_HarvardOxford3mm.nii.gz', 'BG_HarvardOxford3mm.nii.gz', save_path)

"""























