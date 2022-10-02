#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 12:20:39 2022

@author: chvis
"""


import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('default')
import os
import pandas as pd
import numpy as np


## from this website under "calculate morisita Horn" https://github.com/vsarsani/stat-lab/blob/master/MHindex.ipynb

MH_index = []
for x in unique_merged_final.iloc[:,2:]:
    d1 = unique_merged_final[x]
    for y in unique_merged_final.iloc[:,2:]:
        if y.startswith('Library_'):
            d2 = unique_merged_final[y]
            product_d1_d2 = d1*d2
            d1_square = d1*d1
            d2_square = d2*d2
            cal_cols=pd.concat([d1,d2,product_d1_d2,d1_square,d2_square], axis=1).fillna(0)
            cal_cols.columns=['d1','d2','product_d1_d2','d1_square','d2_square']
            sub = cal_cols['product_d1_d2']>0
            subset = cal_cols[sub]
            MH = 2*subset['product_d1_d2'].sum()/(subset['d1_square'].sum()+subset['d2_square'].sum())
            #MH_index.extend([MH])
            MH_index.append(MH)
    
# MH_index_use = MH_index[:]
# var = 1      
# new_df = pd.DataFrame()
# for q in unique_merged_final.iloc[:,2:].columns:
#     new_list = []
#     for val in MH_index_use:
#         print(var)
#         if (var%41)!=0: 
#            new_list.append(val)
#            var=var+1
#         else:
#             new_df[q] = new_list
#             var=var+1
#             break

# MH_index_use = MH_index[:]
 
# new_df = pd.DataFrame()
# for q in unique_merged_final.iloc[:,2:].columns:
#     #print(q)
#     new_list = []
#     var = 0
#     for c in MH_index_use:
#         if (var%41)=0: 
#            new_df[q] = new_list
#            break
#         else: 
#            new_list.append(c)
#            var=var+1

## We tried the below to make this matrix with python, but failed.
## I instead spent 5 minutes to do it manually in excel
## The labels where also added in excel
#%% Dictionaries to easily convert data
## Make a dictionary to put selection rounds into projects
TPL_to_Project = {'TPL0027' : 'Shirin', 
              'TPL0039' : 'Chris_Myo', 
              'TPL0065' : 'Line',
              'TPL0066' : 'Line',
              'TPL0067' : 'Line',
              'TPL0095' : 'Line',
              'TPL0096' : 'Line',
              'TPL0097' : 'Line',
              'TPL0098' : 'Line',
              'TPL0099' : 'Line',
              'TPL0101' : 'Line',
              'TPL0109' : 'Chris',
              'TPL0110' : 'Chris',
              'TPL0111' : 'Chris',
              'TPL0123' : 'Chris',
              'TPL0124' : 'Chris',
              'TPL0125' : 'Chris',
              'TPL0126' : 'Chris',
              'TPL0127' : 'Chris',
              'TPL0130' : 'Chris',
              'TPL0227' : 'Helen',
              'TPL0228' : 'Helen',
              'TPL0229' : 'Helen',
              'TPL0230' : 'Helen',
              'TPL0266' : 'Helen',
              'TPL0204' : 'Isabel',
              'TPL0256' : 'Isabel',
              'TPL0285' : 'Shirin',
              'TPL0298' : 'Shirin?',
              'TPL0301' : 'Shirin',
              'TPL0302' : 'Shirin',
              }



Library_to_TPL = {'Library_1_F01' : 'TPL0049',
             'Library_1_F02' : 'TPL0050',
             'Library_1_F03' : 'TPL0065',
             'Library_1_F05' : 'TPL0066',
             'Library_1_F06' : 'TPL0067',
             'Library_1_F07' : 'TPL0068',
             'Library_1_F08' : 'TPL0095',
             'Library_1_F09' : 'TPL0096',
             'Library_1_F10' : 'TPL0097',
             'Library_1_F11' : 'TPL0098',
             'Library_1_F12' : 'TPL0099',
             'Library_2_F01' : 'TPL0101',
             'Library_2_F02' : 'TPL0103',
             'Library_2_F03' : 'TPL0109',
             'Library_2_F05' : 'TPL0110',
             'Library_2_F06' : 'TPL0111',
             'Library_2_F07' : 'TPL0123',
             'Library_2_F08' : 'TPL0124',
             'Library_2_F09' : 'TPL0125',
             'Library_2_F10' : 'TPL0126',
             'Library_2_F11' : 'TPL0127',
             'Library_2_F12' : 'TPL0128',
             'Library_3_F01' : 'TPL0214',
             'Library_3_F02' : 'TPL0220',
             'Library_3_F03' : 'TPL0221',
             'Library_3_F05' : 'TPL0227',
             'Library_3_F06' : 'TPL0266',
             'Library_3_F07' : 'TPL0228',
             'Library_3_F08' : 'TPL0229',
             'Library_3_F09' : 'TPL0230',
             'Library_3_F10' : 'TPL0037',
             'Library_3_F11' : 'TPL0038',
             'Library_3_F12' : 'TPL0039',
             'Library_4_F01' : 'TPL0129',
             'Library_4_F02' : 'TPL0130',
             'Library_4_F03' : 'TPL0204',
             'Library_4_F05' : 'TPL0239',
             'Library_4_F06' : 'TPL0256',
             'Library_4_F07' : 'TPL0285',
             'Library_4_F08' : 'TPL0301',
             'Library_4_F09' : 'TPL0302',   
    }

Library_to_panning = {'Library_1_F01' : 'E--_K',
             'Library_1_F02' : 'D--_K',
             'Library_1_F03' : 'ED-_K',
             'Library_1_F05' : 'EE-_K',
             'Library_1_F06' : 'DD-_K',
             'Library_1_F07' : 'DE-_K',
             'Library_1_F08' : 'EDE_K',
             'Library_1_F09' : 'EED_K',
             'Library_1_F10' : 'EEE_K',
             'Library_1_F11' : 'DDD_K',
             'Library_1_F12' : 'DDE_K',
             'Library_2_F01' : 'DEE_K',
             'Library_2_F02' : 'A--_K',
             'Library_2_F03' : 'AA-_K',
             'Library_2_F05' : 'AB-_K',
             'Library_2_F06' : 'AC-_K',
             'Library_2_F07' : 'AAA_K',
             'Library_2_F08' : 'AAB_K',
             'Library_2_F09' : 'AAC_K',
             'Library_2_F10' : 'ABA_K',
             'Library_2_F11' : 'ABB_K',
             'Library_2_F12' : 'ABC_K',
             'Library_3_F01' : 'F--_L',
             'Library_3_F02' : 'FF-_L',
             'Library_3_F03' : 'FG-_L',
             'Library_3_F05' : 'FFF_L_100nM',
             'Library_3_F06' : 'FFF_L_10nM',
             'Library_3_F07' : 'FFG_L',
             'Library_3_F08' : 'FGF_L',
             'Library_3_F09' : 'FGG_L',
             'Library_3_F10' : 'B--_L',
             'Library_3_F11' : 'BB-_L',
             'Library_3_F12' : 'BBB_L',
             'Library_4_F01' : 'ACA_K',
             'Library_4_F02' : 'ACB_K',
             'Library_4_F03' : 'H--_L',
             'Library_4_F05' : 'HH-_L',
             'Library_4_F06' : 'HHH_L',
             'Library_4_F07' : 'J--_L_100uM',
             'Library_4_F08' : 'JJ-_L_1nM',
             'Library_4_F09' : 'JJJ_L_200pM',   
    }


Library_2_to_panning = {'Library_1_F01_2' : 'E--_K',
             'Library_1_F02_2' : 'D--_K',
             'Library_1_F03_2' : 'ED-_K',
             'Library_1_F05_2' : 'EE-_K',
             'Library_1_F06_2' : 'DD-_K',
             'Library_1_F07_2' : 'DE-_K',
             'Library_1_F08_2' : 'EDE_K',
             'Library_1_F09_2' : 'EED_K',
             'Library_1_F10_2' : 'EEE_K',
             'Library_1_F11_2' : 'DDD_K',
             'Library_1_F12_2' : 'DDE_K',
             'Library_2_F01_2' : 'DEE_K',
             'Library_2_F02_2' : 'A--_K',
             'Library_2_F03_2' : 'AA-_K',
             'Library_2_F05_2' : 'AB-_K',
             'Library_2_F06_2' : 'AC-_K',
             'Library_2_F07_2' : 'AAA_K',
             'Library_2_F08_2' : 'AAB_K',
             'Library_2_F09_2' : 'AAC_K',
             'Library_2_F10_2' : 'ABA_K',
             'Library_2_F11_2' : 'ABB_K',
             'Library_2_F12_2' : 'ABC_K',
             'Library_3_F01_2' : 'F--_L',
             'Library_3_F02_2' : 'FF-_L',
             'Library_3_F03_2' : 'FG-_L',
             'Library_3_F05_2' : 'FFF_L_100nM',
             'Library_3_F06_2' : 'FFF_L_10nM',
             'Library_3_F07_2' : 'FFG_L',
             'Library_3_F08_2' : 'FGF_L',
             'Library_3_F09_2' : 'FGG_L',
             'Library_3_F10_2' : 'B--_L',
             'Library_3_F11_2' : 'BB-_L',
             'Library_3_F12_2' : 'BBB_L',
             'Library_4_F01_2' : 'ACA_K',
             'Library_4_F02_2' : 'ACB_K',
             'Library_4_F03_2' : 'H--_L',
             'Library_4_F05_2' : 'HH-_L',
             'Library_4_F06_2' : 'HHH_L',
             'Library_4_F07_2' : 'J--_L_100uM',
             'Library_4_F08_2' : 'JJ-_L_1nM',
             'Library_4_F09_2' : 'JJJ_L_200pM',   
    }


TPL_to_panning = {'TPL0049' : 'E--_K',
             'TPL0050' : 'D--_K',
             'TPL0065' : 'ED-_K',
             'TPL0066' : 'EE-_K',
             'TPL0067' : 'DD-_K',
             'TPL0068' : 'DE-_K',
             'TPL0095' : 'EDE_K',
             'TPL0096' : 'EED_K',
             'TPL0097' : 'EEE_K',
             'TPL0098' : 'DDD_K',
             'TPL0099' : 'DDE_K',
             'TPL0101' : 'DEE_K',
             'TPL0103' : 'A--_K',
             'TPL0109' : 'AA-_K',
             'TPL0110' : 'AB-_K',
             'TPL0111' : 'AC-_K',
             'TPL0123' : 'AAA_K',
             'TPL0124' : 'AAB_K',
             'TPL0125' : 'AAC_K',
             'TPL0126' : 'ABA_K',
             'TPL0127' : 'ABB_K',
             'TPL0128' : 'ABC_K',
             'TPL0214' : 'F--_L',
             'TPL0220' : 'FF-_L',
             'TPL0221' : 'FG-_L',
             'TPL0227' : 'FFF_L_100nM',
             'TPL0266' : 'FFF_L_10nM',
             'TPL0228' : 'FFG_L',
             'TPL0229' : 'FGF_L',
             'TPL0230' : 'FGG_L',
             'TPL0037' : 'B--_L',
             'TPL0038' : 'BB-_L',
             'TPL0039' : 'BBB_L',
             'TPL0129' : 'ACA_K',
             'TPL0130' : 'ACA_K',
             'TPL0204' : 'H--_L',
             'TPL0239' : 'HH-_L',
             'TPL0256' : 'HHH_L',
             'TPL0285' : 'J--_L_100uM',
             'TPL0301' : 'JJ-_L_1nM',
             'TPL0302' : 'JJJ_L_200pM',   
    }


#%%
import seaborn as sns
os.chdir('/Users/chvis/Jupyter-Demo/pandas_test')
MH_matrix = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/MH_index_manual.csv', sep=';', index_col=0)
MH_matrix = MH_matrix.drop(index=['Library_2_F10_2', 'Library_3_F08_2', 'Library_3_F09_2']) #L2_F10=19 , L3_F08=28 , L3_09=29
MH_matrix = MH_matrix.drop(columns=['Library_2_F10_2', 'Library_3_F08_2', 'Library_3_F09_2'])
#Color the samples based on libraries and projects
tags = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/Experiments_ex_info_fail_delete.csv',sep=';', index_col=(0))
df_cluster = tags.join(MH_matrix)
df_cluster.rename(columns=Library_2_to_panning, inplace=True) # convert to 
df_cluster.rename(index=Library_2_to_panning, inplace=True)
Lib_color = df_cluster.Lib_ID.map({
    'L1': 'rosybrown',
    'L2': 'indianred',
    'L3': 'maroon',
    'L4': 'red'
})

Project_color = df_cluster.Project.map({
    'Chris': 'slategrey',
    'Line': 'lightsteelblue',
    'Helen': 'cornflowerblue',
    'C_MyoII': 'royalblue',
    'Isabel': 'navy',
    'Shirin': 'blue'
})
#Clustermap used as heatmap by turning off the clustering row_cluster=False, col_cluster=False
numerical_only = df_cluster.columns[2:]
heat_map = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
ax = heat_map.ax_heatmap
ax.plot([11, 11], [0, 741], 'r-', lw = 3) # make a vertical line after L1-L2 seperator
ax.plot([21, 21], [0, 741], 'r-', lw = 3) #L2-L3 seperator
ax.plot([30, 30], [0, 741], 'r-', lw = 3) #L3-L4 seperator
ax.plot([0, 741], [11, 11], 'r-', lw = 3) #L1-L2 seperator 
ax.plot([0, 741], [21, 21], 'r-', lw = 3) #L2-L3 seperator
ax.plot([0, 741], [30, 30], 'r-', lw = 3) #L3-L4 seperator
plt.savefig('heatmap_MH.png', dpi=400)


## clustermap does not really seem relevant
# numerical_only = df_cluster.columns[2:]
# modify = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
# plt.savefig('clustermap_MH.png', dpi=400)
#modify.ax_row_dendrogram.set_visible(False)
### HERE IS THE ISSUE
#modify.ax_heatmap.xaxis.set_label_position('top')
#modify.ax_heatmap.yaxis.set_label_position('left')

           


        

  


            
            


##  Just use for i in range(len(single_row) and step_size_41) and append it to new matrix (list of lists)






### I found this script to make a list into a matrix
### source: https://stackoverflow.com/questions/66049911/converting-list-into-matrix

# L = MH_index
# v = 41

# M = [[] for _ in range(v)]
# r = 0
# for i in L:
#     M[r].append(i)
#     r += 1
#     if r == v:
#         r = 0
#     matrix_MH_index.append(M)


