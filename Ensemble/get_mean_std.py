import nibabel as nib
import numpy as np
import sys
import os


cwd = '/data3/aki/HNU/Refdata/'
subs = os.listdir(cwd)
img = nib.load(cwd+subs[0])
img_data = img.get_fdata()
shape = np.mean(img_data, axis=3).shape

mean_agg = np.zeros([shape[0],shape[1],shape[2],len(subs)])


for s in range(30):
    img = nib.load(cwd+subs[s])
    img_data = img.get_fdata()
    m = np.mean(img_data, axis=3)
    mean_agg[:,:,:,s] = m
    
affine = img.affine
header = img.header


    
mean = np.mean(mean_agg, axis=3)
std = np.std(mean_agg, axis=3)



mean_img = nib.Nifti1Image(mean, affine, header=header)
mean_img.get_data_dtype() == np.dtype(np.int16)

std_img = nib.Nifti1Image(std, affine, header=header)
std_img.get_data_dtype() == np.dtype(np.int16)








