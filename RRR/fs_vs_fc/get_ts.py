# Please run this script with following format:
# python get_ts.py input_folder output_folder time_series_type





from nilearn import image
import nilearn
from nilearn.input_data import NiftiLabelsMasker
from nilearn import datasets
import pickle as pkl


dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
atlas_filename = dataset.maps
labels = dataset.labels

def extract_ts(path, atlasPath):
    img = image.load_img(path)
    template = nilearn.datasets.load_mni152_template()
    img_re = image.resample_to_img(img, template, interpolation='linear')
    masker = NiftiLabelsMasker(labels_img=atlasPath, standardize=True,
                               memory='nilearn_cache', verbose=5)
    ts = masker.fit_transform(img_re)
    return ts


import sys
import os


cwd = sys.argv[1]
files = os.listdir(cwd)

for file in files:
    ts = extract_ts(cwd+file, atlas_filename)
    f = open(sys.argv[2]+file[0:7]+'_' + sys.argv[3]+'.pkl', "wb")
    pkl.dump(ts, f)
    f.close()
