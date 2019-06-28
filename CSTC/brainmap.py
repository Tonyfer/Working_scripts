from nilearn import datasets
from nilearn import image
import nibabel as nib
from nilearn import plotting
import matplotlib.pyplot as plt
import numpy as np



ss = nib.freesurfer.io.read_annot('../documents/rh.aparc.annot')
fsaverage = datasets.fetch_surf_fsaverage()


import os

files = [i for i in [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser('/Users/yu.tong/clusters_analysis/cluster_ana/')) for f in fn] if ('variable.txt' in i) and ('ru' in i)]

for file in files[0]:
    fold = file.replace('variable.txt','')
    f = open(file)
    r_li = [i.strip('\n') for i in f.readlines()]
    f.close()
    f = open(fold+'r_square.txt')
    r_s = ''
    for i in f.readlines():
        r_s += i.replace('\n', ' ')
    f.close()


    v_s_l = list(set([i.split('.')[-1] for i in r_li if '_v_' in i and 'lh' in i]))
    v_s_r = list(set([i.split('.')[-1] for i in r_li if '_v_' in i and 'rh' in i]))
    v_s_l = [v_s_l[i] for i in range(len(v_s_l)) if i%2==0]
    v_s_r = [v_s_r[i] for i in range(len(v_s_r)) if i%2==0]


    regi_l = [i.decode('UTF-8') for i in ss[2]]


    tmp_l = np.zeros(ss[0].shape)
    for i in v_s_l:
        if i in regi_l:
            tmp_l[ss[0]==regi_l.index(i)] = regi_l.index(i)

    tmp_r = np.zeros(ss[0].shape)

    for i in v_s_r:
        if i in regi_l:
            tmp_r[ss[0]==regi_l.index(i)] = regi_l.index(i)
                

    j=1
    plt.figure(figsize=(10, 20))
    for i in ['lateral', 'medial', 'dorsal', 'ventral', 'anterior', 'posterior']:
        ax = plt.subplot(6,2,2*j-1, projection='3d')
        plotting.plot_surf_roi(fsaverage['pial_left'], roi_map=tmp_l,
                           hemi='left', view=i,title='Left_'+i,
                           bg_map=fsaverage['sulc_left'], bg_on_data=True,axes=ax,
                           darkness=1)
        ax = plt.subplot(6,2,2*j, projection='3d')
        plotting.plot_surf_roi(fsaverage['pial_right'], roi_map=tmp_r,
                           hemi='right', view=i,title='Right_'+i,
                           bg_map=fsaverage['sulc_right'], bg_on_data=True,axes=ax,
                           darkness=1)
        j+=1
    plt.savefig(fold+'surface_volume.png')
    
    
    v_sub_l = [i.replace('fs_v_','') for i in r_li if '_v_' in i and ('Left' in i or 'Right' in i)]
    v_sub_l = [i.replace('ant_v_','') for i in v_sub_l]
    v_sub_l = [i.replace('.area','') for i in v_sub_l]
    v_sub_l = list(set([i.replace('.',' ') for i in v_sub_l]))


    harvard_02 = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr0-2mm')
    data = nib.load(harvard_02.maps).get_fdata()
    plot_d = np.zeros(data.shape)
    for i in v_sub_l:
        if i in harvard_02.labels:
            plot_d[data==harvard_02.labels.index(i)] = harvard_02.labels.index(i)
            
    
            
    affine = nib.load(harvard_02.maps).affine
    header = nib.load(harvard_02.maps).header
    plo = nib.Nifti1Image(plot_d, affine, header=header)
    
    plt.figure(figsize=(20, 10))
    if len(np.unique(plot_d))==1:
        template=datasets.load_mni152_template()
        plotting.plot_anat(template, black_bg=False)
        plt.savefig(fold+'sub_volume.png')
    else:
        plotting.plot_roi(plo) 
        plt.savefig(fold+'sub_volume.png')


    
    






for file in files:
    fold = file.replace('variable.txt','')
    f = open(file)
    r_li = [i.strip('\n') for i in f.readlines()]
    f.close()


    a_s_l = list(set([i.split('.')[-1] for i in r_li if 'area' in i and 'lh' in i]))
    a_s_r = list(set([i.split('.')[-1] for i in r_li if 'area' in i and 'rh' in i]))
    a_s_l = [a_s_l[i] for i in range(len(a_s_l)) if i%2==0]
    a_s_r = [a_s_r[i] for i in range(len(a_s_r)) if i%2==0]


    regi_l = [i.decode('UTF-8') for i in ss[2]]


    tmp_l = np.zeros(ss[0].shape)
    for i in a_s_l:
        if i in regi_l:
            tmp_l[ss[0]==regi_l.index(i)] = regi_l.index(i)

    tmp_r = np.zeros(ss[0].shape)

    for i in a_s_r:
        if i in regi_l:
            tmp_r[ss[0]==regi_l.index(i)] = regi_l.index(i)
                

    j=1
    plt.figure(figsize=(10, 20))
    for i in ['lateral', 'medial', 'dorsal', 'ventral', 'anterior', 'posterior']:
        ax = plt.subplot(6,2,2*j-1, projection='3d')
        plotting.plot_surf_roi(fsaverage['pial_left'], roi_map=tmp_l,
                           hemi='left', view=i,title='Left_'+i,
                           bg_map=fsaverage['sulc_left'], bg_on_data=True,axes=ax,
                           darkness=1)
        ax = plt.subplot(6,2,2*j, projection='3d')
        plotting.plot_surf_roi(fsaverage['pial_right'], roi_map=tmp_r,
                           hemi='right', view=i,title='Right_'+i,
                           bg_map=fsaverage['sulc_right'], bg_on_data=True,axes=ax,
                           darkness=1)
        j+=1
    plt.savefig(fold+'surface_volume.png')






    
    
    









