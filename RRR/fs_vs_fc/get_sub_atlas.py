from nilearn import datasets
from nilearn import image
ho = datasets.fetch_atlas_harvard_oxford('sub-prob-2mm', data_dir=None, symmetric_split=False, resume=True, verbose=1)

lables = ['Left Caudate', 'Right Caudate', 'Left Putamen', 'Right Putamen', 'Left Thalamus', 'Right Thalamus', 'Left Accumbens', 'Right Accumbens', 'Left Pallidum', 'Right Pallidum']

indexs = [ho.labels.index(i)-1 for i in lables]

import nibabel as nib
img = nib.load(ho.maps)

ho_prob = img.get_fdata()

affine = img.affine
header = img.header


import numpy as np
at = np.zeros([91,109,91])

for i in indexs:
    tt = ho_prob[:,:,:,i]
    at[tt>33] = i+1
    

at33 = nib.Nifti1Image(at, affine, header=header)
at33.get_data_dtype() == np.dtype(np.int16)

nib.nifti1.save(at33, 'subcortical_33.nii.gz')




















