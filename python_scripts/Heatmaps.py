#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 11:08:50 2022

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
#%% RUN ONLY WHEN NEED TO OBTAIN DATA - Getting the data: Matching values between datasets - 2 hour run time - nucleotides

# What i need to do is first to construct a pairwise overlap matrix where each unique clone is 
# from each dataset is compared to all clones from all other datasets - When there is a match, 
# this gives a score of 1. Then all matches between samples are counted and displayed as a 
# pairwise overlap matrix with heatmaps - on one side of the matrix use nt-seqs, and AA-seqs 
# on the other side. - This analysis just not take cloneCounts/fraction into consideration - 
# which means that anything with a cloneCount above 1 will be equal to each other.

# Next we want to get it quantified using Morisita Horn, an example for preparing the data for 
# this can be seen here

# We have to prepare the data for morisita horn
# example
# Set1: Car, CTR, CSR
# Set2: Car, CLR, CWR
# Here the overlap would be 1/3
# Using Morisita Horn
# Set1 - Set2
#  CAR 50% - CAR 40%
#  CTR 30% - CLR 40%
#  CSR 20% - CWR 20%
# Rank them by the CDRs
#  CAR 50% - 40%
#  CTR 30% - 0%
#  CSR 20% - 0
#  CLR 0% - 40% CWR 0% - 20%
#  To get this in R you merge
# This matrix you give to Morisita Horn
# Potential fixes

# len(df)-len(df.drop_duplicates()) for a quick fix to counting identical overlaps between samples. 
# Won't help for further processing i believe. I think the groupby function could make the preparation 
# for Morisita Horn easier


## This script takes 2 hours to run, therefore it has been commented out
## When you want to run it, remove the # from the four lines below here
## This script extract the nucleotide sequences and counts how many matches are between them and different samples
finaldataframe = pd.DataFrame()
addedlibcolumn = False #Checks whether the "Libs" column has been filled
libsnamestoadd = [] #Used to add the libraries to the "Libs" column
for Match_from_list in list1:
    arraydata = []
    dataframedata = []
    if addedlibcolumn == False:
        finaldataframe['Libs'] = libsnamestoadd

    if Match_from_list.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
        Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+Match_from_list, sep='\t')
        Export = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
        Export = Export[(Export['cloneCount'] > 1)]   # by running this i get rid of the clones with only 1 cloneCount
        nt_sequences = (Export.nSeqCDR3)
        df_Standard = pd.DataFrame(nt_sequences) #Make into actual dataframe - can be checked with ifinstance funciton
        
        for Match_to_list in list1:
                
            if Match_to_list.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
                if addedlibcolumn == False:
                     libsnamestoadd.append(Match_to_list[0:Match_to_list.index('_UniqueCDR3_Exp_UniqueCDR3.txt')])
                print([Match_to_list[0:Match_to_list.index('_UniqueCDR3_Exp_UniqueCDR3.txt')]])
                Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+Match_to_list, sep='\t')
                Export = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa transeslation)
                Export = Export[(Export['cloneCount'] > 1)]   # by running this i get rid of the clones with only 1 cloneCount
                nt_sequences = (Export.nSeqCDR3) # here aminoacids!!!
                df_nt_sequences = pd.DataFrame(nt_sequences) #Make into actual dataframe - can be checked with ifinstance funciton
                for col in df_nt_sequences:
                    nt_matches = 0
                    for j, row_value in df_nt_sequences[col].iteritems():
                        found = len(df_Standard[df_Standard[col] == df_nt_sequences[col][j]])
                        nt_matches = nt_matches + found
                    #arraydata.append([Match_to_list[0:Match_to_list.index('_UniqueCDR3_Exp_UniqueCDR3.txt')], nt_matches])
                #if Match_to_list != Match_from_list:
                    dataframedata.append({nt_matches})
                print(dataframedata)
        if addedlibcolumn == False:
            finaldataframe['Libs'] = libsnamestoadd # Filles the empty 'Libs'column with the library names
            addedlibcolumn = True
        
        tempdataframe = pd.DataFrame(dataframedata)
        finaldataframe[Match_from_list[0:Match_from_list.index('_UniqueCDR3_Exp_UniqueCDR3.txt')]] = tempdataframe
        print(finaldataframe)
        # finaldataframe = pd.merge(finaldataframe, tempdataframe, on=Match_from_list)
        
#finaldataframe.to_excel("finaldataframe_nt.xlsx")
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
    sample_name.index=sample_name.index.str.replace('_F1_', '_F01_')
    sample_name.index=sample_name.index.str.replace('_F2_', '_F02_')
    sample_name.index=sample_name.index.str.replace('_F3_', '_F03_')
    sample_name.index=sample_name.index.str.replace('_F5_', '_F05_')
    sample_name.index=sample_name.index.str.replace('_F6_', '_F06_')
    sample_name.index=sample_name.index.str.replace('_F7_', '_F07_')
    sample_name.index=sample_name.index.str.replace('_F8_', '_F08_')
    sample_name.index=sample_name.index.str.replace('_F9_', '_F09_')
#%% converting and sorting full nucleotide dataset - post remove 1's and nt seqs not divisible by 3
os.chdir('/Users/chvis/Jupyter-Demo/pandas_test') 
Nucleotide_match = pd.read_excel('finaldataframe_nt.xlsx',index_col=(1)).drop(labels="Unnamed: 0", axis=1)
convert_F1_F01_etc(Nucleotide_match)
Nucleotide_match = Nucleotide_match.sort_index(axis=0)
Nucleotide_match = Nucleotide_match.sort_index(axis=1)
Nucleotide_match = Nucleotide_match.drop(index=['Library_2_F10_2','Library_3_F08_2', 'Library_3_F09_2'], columns=['Library_2_F10_2','Library_3_F08_2', 'Library_3_F09_2']) #drop the failed sequences

Nucleotide_match.to_csv('Nucleotide_match.csv')

#%% converting and sorting top100

finaldataframe_top100_sorted = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/finaldataframe_top100.csv',index_col=(1)).drop(labels="Unnamed: 0", axis=1)
convert_F1_F01_etc(finaldataframe_top100_sorted)
finaldataframe_top100_sorted = finaldataframe_top100_sorted.sort_index(axis=0)
finaldataframe_top100_sorted = finaldataframe_top100_sorted.sort_index(axis=1)
finaldataframe_top100_sorted = finaldataframe_top100_sorted.drop(index=['Library_2_F10_2','Library_3_F08_2', 'Library_3_F09_2'], columns=['Library_2_F10_2','Library_3_F08_2', 'Library_3_F09_2']) #drop the failed sequences

finaldataframe_top100_sorted.to_csv('finaldataframe_top100_sorted.csv')   # only uncomment this when you want to export


tags = pd.read_excel('/Users/chvis/Jupyter-Demo/pandas_test/Experiments_extra_information.xlsx',index_col=(0))

#join the matrix with the file containing library labels and project labels 
top100_sort_tags = tags.join(finaldataframe_top100_sorted)
#top100_sort_tags.to_csv('top100_sort_tags.csv')
#finaldataframe_top100_sorted.to_csv('finaldataframe_top100_sorted.csv')
#%% HEATMAP_ALL       -      get puneet to explain the "for i inxxx and for j in xxx"

os.chdir('/Users/chvis/Jupyter-Demo/pandas_test') 
df_heat = pd.read_csv('df_combined3.csv').drop(labels="Unnamed: 0", axis=1)
df_heat = df_heat.drop(index=[19, 28, 29]) #L2_F10=19 , L3_F08=28 , L3_09=29
df2_heat = pd.read_csv('Nucleotide_match.csv',index_col=(0))
df3_heat = np.zeros([38,38]) #change to dataset size
for i in range(len(df3_heat)):
    for j in range(len(df3_heat)):
        df3_heat[i,j]= df2_heat.iloc[i,j]/((df_heat.loc[df_heat['Experiment']==df2_heat.index[i], "Uniqueclones_>1"].iloc[0]+df_heat.loc[df_heat['Experiment']==df2_heat.columns[j], "Uniqueclones_>1"].iloc[0])/2)

df4_heat = pd.DataFrame(df3_heat, columns= df2_heat.columns, index= df2_heat.index)

#Color the samples based on libraries and projects
tags = pd.read_excel('/Users/chvis/Jupyter-Demo/pandas_test/Experiments_ex_info_fail_delete.xlsx',index_col=(0))
df_cluster = tags.join(df4_heat)


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
heat_map = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
plt.savefig('heatmap_all.png', dpi=400)

numerical_only = df_cluster.columns[2:]
modify = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
plt.savefig('clustermap_all.png', dpi=400)
#modify.ax_row_dendrogram.set_visible(False)
### HERE IS THE ISSUE
#modify.ax_heatmap.xaxis.set_label_position('top')
#modify.ax_heatmap.yaxis.set_label_position('left')
###
#%% HEATMAP_top100       -      get puneet to explain the "for i inxxx and for j in xxx"

os.chdir('/Users/chvis/Jupyter-Demo/pandas_test') 
df_heat = pd.read_csv('df_combined3.csv').drop(labels="Unnamed: 0", axis=1)
df_heat = df_heat.drop(index=[19, 28, 29]) #L2_F10=19 , L3_F08=28 , L3_09=29
df2_heat = pd.read_csv('finaldataframe_top100_sorted.csv',index_col=(0))
df3_heat = np.zeros([38,38]) #change to dataset size
for i in range(len(df3_heat)):
    for j in range(len(df3_heat)):
        df3_heat[i,j]= df2_heat.iloc[i,j]/((df_heat.loc[df_heat['Experiment']==df2_heat.index[i], "Uniqueclones_>1"].iloc[0]+df_heat.loc[df_heat['Experiment']==df2_heat.columns[j], "Uniqueclones_>1"].iloc[0])/2)


df4_heat = pd.DataFrame(df3_heat, columns= df2_heat.columns, index= df2_heat.index)

#Color the samples based on libraries and projects
tags = pd.read_excel('/Users/chvis/Jupyter-Demo/pandas_test/Experiments_ex_info_fail_delete.xlsx',index_col=(0))
df_cluster = tags.join(df4_heat)


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
heat_map = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
plt.savefig('heatmap_top100.png', dpi=400)

#numerical_only = df_cluster.columns[2:]
cluster_map = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
plt.savefig('clustermap_top100.png', dpi=400)
#modify.ax_row_dendrogram.set_visible(False)  cluster still on, but dendrogram removed
### HERE IS THE ISSUE
#modify.ax_heatmap.xaxis.set_label_position('top')
#modify.ax_heatmap.yaxis.set_label_position('left')
#%% TEST WITH RAW top100 values
df4_heat = pd.read_csv('finaldataframe_top100_sorted.csv',index_col=(0))

#Color the samples based on libraries and projects
tags = pd.read_excel('/Users/chvis/Jupyter-Demo/pandas_test/Experiments_ex_info_fail_delete.xlsx',index_col=(0))
df_cluster = tags.join(df4_heat)


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
heat_map = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", row_cluster=False, col_cluster=False, annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
plt.savefig('heatmap_top100_raw.png', dpi=400)

#numerical_only = df_cluster.columns[2:]
cluster_map = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), cmap="rocket", annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
plt.savefig('clustermap_top100_raw.png', dpi=400)