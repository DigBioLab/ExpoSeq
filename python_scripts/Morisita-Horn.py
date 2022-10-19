#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 16:15:00 2022

@author: chvis
"""
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('default')
import numpy as np
import fnmatch as fn
import seaborn as sns
import os
#%% RUN EVERTIME - Define converting tool 

##converting F1 -> F01 etc - not sure which way is best
##This is how to add a zero to before the single digit Fs (F1, F2 etc) But it does not work in column rows..
## Should work in both columns and indexes by using df.column.replace / df.index.replace
#Two_digit_back_up_df = back_up_df.replace(regex=['F1_2', 'F2_2', 'F3_2', 'F5_2', 'F6_2', 'F7_2', 'F8_2', 'F9_2'], value=['F01_2', 'F02_2', 'F03_2', 'F05_2', 'F06_2', 'F07_2', 'F08_2', 'F09_2'])

def convert_F1_F01_etc(sample_name): # by typing in sample name i.e. unique_merged - all columns with _F1_ will be convertedd to _F01_ etc
    sample_name.columns=sample_name.columns.str.replace('_F1_', '_F01_')
    sample_name.columns=sample_name.columns.str.replace('_F2_', '_F02_')
    sample_name.columns=sample_name.columns.str.replace('_F3_', '_F03_')
    sample_name.columns=sample_name.columns.str.replace('_F5_', '_F05_')
    sample_name.columns=sample_name.columns.str.replace('_F6_', '_F06_')
    sample_name.columns=sample_name.columns.str.replace('_F7_', '_F07_')
    sample_name.columns=sample_name.columns.str.replace('_F8_', '_F08_')
    sample_name.columns=sample_name.columns.str.replace('_F9_', '_F09_')

#%% Morisita-Horn - make the combined list of all unique sequences in all dataset with 
all_seq_combined = []
for allseq in list1:
    if allseq.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
        Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+allseq, sep='\t')
        Export = Export[(Export['cloneCount'] > 1)] 
        Export = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
        nt_stack_input = (Export.loc[0:10,'nSeqCDR3']).to_list()
        all_seq_combined.extend(nt_stack_input)
        df_stacked = pd.DataFrame(nt_stack_input)
seq = []
freq =[]
for column in df_stacked:
    [seq.append(x) for x in df_stacked [column] if x not in seq] #append sequences to seq if not found in seq already. How about a partial match?

numpy_list = np.array(all_seq_combined)
(unique,count) = np.unique(numpy_list, return_counts=True)
frequncies = np.asarray((unique,count)).T
unique_nt_df = pd.DataFrame(frequncies, columns = ['nSeqCDR3','count']) # count?
#print(df)   
#unique_nt_df['count'].hist(bins=100)#column =['CDR','count'])
print(unique_nt_df)

#unique_nt_df.to_excel('unique_nt_df.xlsx')

## Merging them
unique_merged = unique_nt_df[:]
for match_to in list1:
    if match_to.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
        Export = pd.read_csv(match_to, sep='\t')
        Export1 = Export[(Export['cloneCount'] > 1)] 
        Export2 = Export1[(Export1['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)

        column = Export2.groupby('nSeqCDR3', as_index=False)['cloneFraction'].cumsum()
        Export2['cloneFraction'] = column['cloneFraction']
        Export2 = Export2.drop_duplicates(subset= 'nSeqCDR3')

        match_input = (Export2.loc[0:10,['nSeqCDR3', 'cloneFraction']])
        unique_merged = unique_merged.merge(match_input,how='left', on='nSeqCDR3')
        unique_merged = unique_merged.rename(columns={'cloneFraction': match_to[0:match_to.index('_UniqueCDR3_Exp_UniqueCDR3.txt')]})
print(unique_merged)
convert_F1_F01_etc(unique_merged)
unique_merged=unique_merged.sort_index(axis=1)



