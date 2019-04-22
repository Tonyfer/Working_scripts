import pandas as pd
import numpy as np
import sys
import os

name_v_ant_template=[]
with open('/data3/cdb/ytong/merged_tables/templates/v_ant.txt','r') as f:
    for line in f:
        name_v_ant_template+=(line.strip('\n').split(','))

        
name_v_fs_template=[]

with open('/data3/cdb/ytong/merged_tables/templates/v_fs.txt','r') as f:
	for line in f:
		name_v_fs_template+=(line.strip('\n').split(','))
        
        

        

def get_df_cbic(cwd, tech, group):
    def get_v_a(s):
        if s in name_v_ant:
            return(volum_ant[volum_ant.columns[2]][name_v_ant.index(s)])
        else:
            return ''

    def get_v_fs(s):
        if s in name_v_fs:
            return(volum_fs[volum_fs.columns[2]][name_v_fs.index(s)])
        else:
            return ''
        
        
    files = os.listdir(cwd)
    fi_list = []
    for name in files:
        if 'sub' in name:
            fi_list.append(name)
    f_list = []
    for path in fi_list:
        a = os.path.exists(cwd+path+'/'+tech+'/mindboggle/')
        b = os.path.exists(cwd+path+'/'+tech+'/mindboggle/tables/left_cortical_surface/')
        c = os.path.exists(cwd+path+'/'+tech+'/mindboggle/tables/right_cortical_surface/')
        f = os.path.exists(cwd+path+'/'+tech+'/mindboggle/tables/volume_per_freesurfer_label.csv')
        g = os.path.exists(cwd+path+'/'+tech+'/mindboggle/tables/volume_per_ants_label.csv')
        if a and b and c and f and g:
            f_list.append(path)
    lh = pd.read_csv(cwd+f_list[0]+'/'+tech+'/mindboggle/tables/left_cortical_surface/label_shapes.csv')
    name_s_a = ['area' + lh['name'][i] for i in range(len(lh['name']))]
    name_s_tdm = ['tdm' + lh['name'][i] for i in range(len(lh['name']))]
    name_s_t = ['thickness'+lh['name'][i] for i in range(len(lh['name']))]
    #name_t_fs = ['fs_t_'+i for i in name_t_fs_template]
    name_v_fs = ['fs_v_'+i for i in name_v_fs_template]
    #name_t_ant = ['ant_t_'+i for i in name_t_ant_template]
    name_v_ant = ['ant_v_'+i for i in name_v_ant_template]
    var = name_v_fs+name_v_ant+name_s_a+name_s_tdm+name_s_t
    df = pd.DataFrame(np.zeros([1,len(var)]), columns=var)
    i = 0
    for sub in f_list:
        sub = cwd + sub
        volum_fs = pd.read_csv(sub+'/'+tech+'/mindboggle/tables/volume_per_freesurfer_label.csv')
        volum_ant = pd.read_csv(sub+'/'+tech+'/mindboggle/tables/volume_per_ants_label.csv')
        lh = pd.read_csv(sub+'/'+tech+'/mindboggle/tables/left_cortical_surface/label_shapes.csv')
        rh = pd.read_csv(sub+'/'+tech+'/mindboggle/tables/right_cortical_surface/label_shapes.csv')
        name_v_ant = [volum_ant['name'][i] for i in range(len(volum_ant['name']))]
        value_v_ant = [get_v_a(i) for i in name_v_ant_template]
        #name_v_ant = ['ant_v_'+volum_ant['name'][i] for i in range(len(volum_ant['name']))]
        name_v_fs = [volum_fs['name'][i] for i in range(len(volum_fs['name']))]
        value_v_fs = [get_v_fs(i) for i in name_v_fs_template]
        

        va_l_a = [lh['area'][i] for i in range(len(lh['name']))]
        va_r_a = [rh['area'][i] for i in range(len(rh['name']))]
        #name_s_tdm = ['tdm' + lh['name'][i] for i in range(len(lh['name']))]
        va_l_tdm = [lh['travel depth: median'][i] for i in range(len(lh['name']))]
        va_r_tdm = [rh['travel depth: median'][i] for i in range(len(rh['name']))]
        va_l_th = [lh['freesurfer thickness: median'][i] for i in range(len(lh['name']))]
        va_r_th = [rh['freesurfer thickness: median'][i] for i in range(len(rh['name']))]
        val_a = [va_l_a[i]+va_r_a[i] for i in range(len(va_l_a))]
        val_tdm = [va_l_tdm[i]+va_r_tdm[i] for i in range(len(va_l_tdm))]
        val_th = [va_l_th[i]+va_r_th[i] for i in range(len(va_l_th))]
    
        val = value_v_fs+value_v_ant+val_a+val_tdm+val_th
        df.loc[i] = val
        i+=1
    df['sub'] = f_list
    df['group'] = [group for i in range(len(f_list))]
    df['tech'] = [tech for i in range(len(f_list))]
    return(df)

cbic_list = ['/data2/HBNcore/CMI_HBN_Data/MRI/CBIC/Release_Derivatives/R3_20170801_20171131/data_release_ExternalID/', '/data2/HBNcore/CMI_HBN_Data/MRI/CBIC/Release_Derivatives/R4_20171201_20180531/data_release_ExternalID/', '/data2/HBNcore/CMI_HBN_Data/MRI/CBIC/Release_Derivatives/R5_20180601_20180931/data_release_ExternalID/', '/data2/HBNcore/CMI_HBN_Data/MRI/CBIC/Release_Derivatives/R6_20181001_20181231/data_release_ExternalID/']


tech_list = ['T1w_HCP', 'T1w_VNav', 'T1w_VNavNorm']


data = [get_df_cbic(cbic_list[j], tech_list[i], i) for i in range(3) for j in range(len(cbic_list))]



df_all = pd.concat([data[i] for i in range(len(data))])


df_all.to_csv('cbic_123.csv')








