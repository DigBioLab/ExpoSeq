#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 16:50:20 2022

@author: chvis
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


#%%
def convert_F1_F01_etc(sample_name): # by typing in sample name i.e. unique_merged - all columns with _F1_ will be convertedd to _F01_ etc
    sample_name.columns=sample_name.columns.str.replace('_F1_', '_F01_')
    sample_name.columns=sample_name.columns.str.replace('_F2_', '_F02_')
    sample_name.columns=sample_name.columns.str.replace('_F3_', '_F03_')
    sample_name.columns=sample_name.columns.str.replace('_F5_', '_F05_')
    sample_name.columns=sample_name.columns.str.replace('_F6_', '_F06_')
    sample_name.columns=sample_name.columns.str.replace('_F7_', '_F07_')
    sample_name.columns=sample_name.columns.str.replace('_F8_', '_F08_')
    sample_name.columns=sample_name.columns.str.replace('_F9_', '_F09_')

#%%
list1 = os.listdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/')


input_file = pd.read_csv('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/3_Nt-sequences_scFv.txt', sep='\t')
input_file = input_file[input_file['V-DOMAIN ID'].str.endswith('_H')]
CDR3_only = input_file[['V-DOMAIN ID', 'JUNCTION']]
CDR3_only= CDR3_only.rename(columns={'JUNCTION':'nSeqCDR3'})

#Sang_NGS_merge = CDR3_only.merge(seqs_counts, how='left', on='nSeqCDR3')

## need to figure out how to get library IDs next to the NGS seqs. These will then be match against Sanger
## seqs, so it will be A: NGS library, B: shared sequence, C: Sanger_ID
CDR3_only = CDR3_only.rename(columns={'V-DOMAIN ID': 'TPL_ID'})

## Grooiming the TPL_IDs to an identical fashion
CDR3_only.replace(to_replace='_TPL', value='TPL', regex=True, inplace=True)
CDR3_only['TPL_ID'] = [x[:14] for x in CDR3_only['TPL_ID']] ## (only take the 14 first characters of each ID, thereby we remove the primer info)
## written by Punnet - this corrects the A1 top A01 in the TPL0039 samples. If the last character is _ then insert 0 inbetween A and 1 etc -> A01
CDR3_only['TPL_ID'] = [str(x[:-2]+'0'+x[-2]) if x[-1]=='_' else x for x in CDR3_only['TPL_ID']]

## When matching the sequences - to accomodate for imgt and mixcr including different things, we had to do the following
## mixcr = CARxxxxxxWG
## IMGT = Starts from A and ends before WG
## Therefore to be able to merge on direct matches, i have to remove the first 3 and the last 6 nucleotides
## of the mixcr analysis data

## First i need to find a way to bring all NGS sequences in with their respective library 
## Could i use the matrix setup used for pre-MH analysis
## Instead of allseqs on the left, we will have SANGER seqs, and then each library as already is

CDR3_match = CDR3_only.loc[:,['nSeqCDR3','TPL_ID']]
for match_to in list1:
    if match_to.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
        Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+match_to, sep='\t')
        Export1 = Export[(Export['cloneCount'] > 1)] 
        Export2 = Export1[(Export1['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
        match_input = (Export2.loc[:,['nSeqCDR3', 'cloneFraction']])
        CDR3_match['nSeqCDR3'] = CDR3_match['nSeqCDR3'].str.upper()
        CDR3_match = CDR3_match.merge(match_input,how='left', on='nSeqCDR3')
        CDR3_match = CDR3_match.rename(columns={'cloneFraction': match_to[0:match_to.index('_UniqueCDR3_Exp_UniqueCDR3.txt')]})
print(CDR3_match)
convert_F1_F01_etc(CDR3_match)
CDR3_label = CDR3_match.loc[:,['nSeqCDR3', 'TPL_ID']]
CDR3_Lib_Sort = CDR3_match.iloc[:,2:].sort_index(axis=1)
CDR3_sang_match_final = pd.concat([CDR3_label,CDR3_Lib_Sort], axis=1)
CDR3_sang_match_final = CDR3_sang_match_final.drop(columns=['Library_2_F10_2', 'Library_3_F08_2', 'Library_3_F09_2'])
CDR3_sang_match_final = CDR3_sang_match_final.drop_duplicates(subset = ['TPL_ID']) # remove duplicate TPL_IDs, some of Helens were sequenced twice
CDR3_sang_match_final = CDR3_sang_match_final.dropna(subset=['nSeqCDR3'])
CDR3_sang_match_final = CDR3_sang_match_final.sort_values('TPL_ID')


#CDR3_sang_match_final_num
#heat_map = sns.clustermap(CDR3_sang_match_final.loc[:,'Library_1_F01_2':], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False)
#plt.savefig('Sang_NGS_match.png', dpi=400)

## not useable anymore, but handy nline to have:
## Export2['nSeqCDR3'] = Export2['nSeqCDR3'].str[3:-3] ## reformat mixcr to imgt layout - remove first 3 and last 6

## Script for dividing TPLxxxx into projects
##


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

## Make a dictionary to put Library names into panning names
import pandas as pd

CDR3_sang_match_final.rename(columns=Library_2_to_panning, inplace=True)  # convert using the Library_2_to_panning dictionary

TPL_Sang_Label = pd.DataFrame()
TPL_Sang_Label['TPL_ID'] = CDR3_sang_match_final['TPL_ID']
TPL_Sang_Label['TPL_ID_short'] = [x[:7] for x in CDR3_sang_match_final['TPL_ID']] # put in your dataframe

TPL_Sang_Label['Project'] = TPL_Sang_Label["TPL_ID_short"].map(TPL_to_Project)

## For the data sorted on the libraries
Sang_NGS_heat_libsort = TPL_Sang_Label.merge(CDR3_sang_match_final, on='TPL_ID', )
Sang_NGS_heat_libsort = Sang_NGS_heat_libsort.sort_values('Project')
Sang_NGS_heat_libsort.rename(columns=Library_2_to_panning, inplace=True)

## For the data to be sorted on projects and relabelled with panning names
# Sang_NGS_heat_pan = TPL_Sang_Label.merge(CDR3_sang_match_final, on='TPL_ID', )
# Sang_NGS_heat_pan.rename(columns=Library_2_to_panning, inplace=True)  # convert using the Library_2_to_panning dictionary

# Sang_NGS_heat_pan_1 = Sang_NGS_heat_pan.iloc[:,:4] # split data - this is not to be sorted
# Sang_NGS_heat_pan_2 = Sang_NGS_heat_pan.iloc[:,4:] # split data - this is to be sorted
# Sang_NGS_heat_pan_2 = Sang_NGS_heat_pan_2.sort_index(axis=1) # sort them
# Sang_NGS_heat_pansort = Sang_NGS_heat_pan_1.join(Sang_NGS_heat_pan_2) # join them back together
# Sang_NGS_heat_pansort = Sang_NGS_heat_pansort.sort_values('Project') # and sort on projects

# ## make heatmap of lib sorted

# Project_color = Sang_NGS_heat_libsort.Project.map({
#     'Chris': 'slategrey',
#     'Line': 'lightsteelblue',
#     'Helen': 'cornflowerblue',
#     'Chris_Myo': 'royalblue',
#     'Isabel': 'navy',
#     'Shirin': 'blue',
#     'Shirin?': 'blue'
#})
#heat_map = sns.clustermap(Sang_NGS_heat_libsort.loc[:,'Library_1_F01_2':], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, row_colors=(Project_color))
#plt.savefig('Sang_NGS_heat_libsort.png', dpi=400)

# ## make heatmap of pan sorted
# Project_color = Sang_NGS_heat_pansort.Project.map({
#     'Chris': 'slategrey',
#     'Line': 'lightsteelblue',
#     'Helen': 'cornflowerblue',
#     'Chris_Myo': 'royalblue',
#     'Isabel': 'navy',
#     'Shirin': 'blue',
#     'Shirin?': 'blue'
# })
# heat_map = sns.clustermap(Sang_NGS_heat_pansort.loc[:,'A--_K':], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, row_colors=(Project_color))
# plt.savefig('Sang_NGS_heat_pansort.png', dpi=400)

#%% trying to change the heatmap so numbers are panning specific, thereby each panning will have a total of 1.0 color to give out to its row


## verically normalized graph libsort

df_heat = Sang_NGS_heat_libsort.iloc[:,4:]
temp = pd.DataFrame()
for i in df_heat:
    temp[i] = df_heat[i]/df_heat[i].sum() # convert values to fractions of the total vertical values to get an even numbering of the values (0 to 1)
    print(temp[i])
temp_labels = Sang_NGS_heat_libsort.iloc[:,:4]
Sang_NGS_heat_libsort_vertnorm = temp_labels.join(temp)
#%% I tried making a function from it but it only returned the temp as a print and not as a useable fashio, can be easily fixed
    ### i tried making it into a function, but i won't return temp as a variable in the variable explorer but only print it
# def heatmap_convert_to_vertfractions(df_of_interest): # put in a matrix containing values only and this will convert all values to the value divided by the sum of the column
#     df_heat = df_of_interest
#     temp = pd.DataFrame()
#     for i in df_heat:
#         temp[i] = df_heat[i]/df_heat[i].sum() # convert values to fractions of the total vertical values to get an even numbering of the values (0 to 1)
#     return  pd.DataFrame(temp)             
#%% 

## Try making a graph as above from it.
    
Project_color = Sang_NGS_heat_libsort_vertnorm.Project.map({
    'Chris': 'slategrey',
    'Line': 'lightsteelblue',
    'Helen': 'cornflowerblue',
    'Chris_Myo': 'royalblue',
    'Isabel': 'navy',
    'Shirin': 'blue',
    'Shirin?': 'blue'
})

heat_map = sns.clustermap(Sang_NGS_heat_libsort_vertnorm.loc[:,'E--_K':], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, row_colors=(Project_color), vmax=.05)
# i would prefer to draw lines for each library project etc, but that is not possible so far
ax.plot([11, 11], [0, 741], 'r-', lw = 3) # make a vertical line after sample 11 ranging for all 741 sanger sequences - color = red
ax.plot([21, 21], [0, 741], 'r-', lw = 3)
ax.plot([30, 30], [0, 741], 'r-', lw = 3)
plt.savefig('Sang_NGS_heat_libsort_vertnorm_max05.png', dpi=400)






## verically normalized graph pansort

# df_heat = Sang_NGS_heat_pansort.iloc[:,4:]
# temp = pd.DataFrame()
# for i in df_heat:
#     temp[i] = df_heat[i]/df_heat[i].sum() # convert values to fractions of the total vertical values to get an even numbering of the values (0 to 1)
# temp_labels = Sang_NGS_heat_pansort.iloc[:,:4]
# Sang_NGS_heat_pansort_vertnorm = temp_labels.join(temp)
# #%% I tried making a function from it but it only returned the temp as a print and not as a useable fashio, can be easily fixed
#     ### i tried making it into a function, but i won't return temp as a variable in the variable explorer but only print it
# # def heatmap_convert_to_vertfractions(df_of_interest): # put in a matrix containing values only and this will convert all values to the value divided by the sum of the column
# #     df_heat = df_of_interest
# #     temp = pd.DataFrame()
# #     for i in df_heat:
# #         temp[i] = df_heat[i]/df_heat[i].sum() # convert values to fractions of the total vertical values to get an even numbering of the values (0 to 1)
# #     return  pd.DataFrame(temp)             

# Project_color = Sang_NGS_heat_pansort_vertnorm.Project.map({
#     'Chris': 'slategrey',
#     'Line': 'lightsteelblue',
#     'Helen': 'cornflowerblue',
#     'Chris_Myo': 'royalblue',
#     'Isabel': 'navy',
#     'Shirin': 'blue',
#     'Shirin?': 'blue'
# })

# heat_map = sns.clustermap(Sang_NGS_heat_pansort_vertnorm.loc[:,'A--_K':], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, row_colors=(Project_color), vmax=.05)
# # i would prefer to draw lines for each library project etc, but that is not possible so far
# ax = heat_map.ax_heatmap
# ax.plot([11, 11], [0, 741], 'r-', lw = 3)
# ax.plot([14, 14], [0, 741], 'r-', lw = 3)
# ax.plot([26, 26], [0, 741], 'r-', lw = 3)
# ax.plot([32, 32], [0, 741], 'r-', lw = 3)
# ax.plot([35, 35], [0, 741], 'r-', lw = 3)
# plt.savefig('Sang_NGS_heat_pansort_vertnorm_max05.png', dpi=400)
