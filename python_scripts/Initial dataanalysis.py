#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 18:36:59 2021

@author: chvis
"""
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('default')
import numpy as np
import fnmatch as fn
#import seaborn as sns
import os
#%% ESSENTIAL RUN EVERYTIME - df_combined3 + list1 defined + seqslost and unique clones lost from removing 1-reads

# This section defines df_combined3 and should therefore always be run
# Further, plots such as seqslost and unique seqs lost form deleting 1-reads will are constructed here

## Example on how to grab data
## df = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/Library_1_F3_2.fastq_AlignmentReport.txt',sep=':')

## Set the working directory
list1 = os.listdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/')

#for each element in list - if ends with txt -. give name
[x for x in list1 if x.endswith('txt')] #maybe not necessary


#I am unsure what i need to use the below line for, therefore it is commented out
#[x[0:x.index('_UniqueCDR3_Exp_UniqueCDR3.txt')]for x in list1 if x.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt')]


# # This python_scripts has been changed to the one below, but i am keeping it in case i need to run something similar
# for i in list1:
#     if i.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
#         #print (i+',')
#         Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+i, sep='\t')
#         evalues = (str(Export.cloneCount.sum())+'\n') #\n = go tothe next line
# #print (evalues)


## For mixcr output do the following
example_file = "C://Users//nils//Desktop//Tropical Pharmaceutical Lab//Project//Library_1_F1_2_UniqueCDR3_Exp_UniqueCDR3.txt"
output = []
for i in list1:
    if i.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'): ## pick the ending that is on your CDR3 files
        Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+i, sep='\t') ## Import all seqs
        Export2 = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
        evalues = (str(Export2.cloneCount.sum())) #\n = go to the next line
        output.append([i[0:i.index('_UniqueCDR3_Exp_UniqueCDR3.txt')], evalues])
#print (output)




# find identical amino acid strangs and add
# include in above loop
# delete the copies
# add count and fraction to the amino acid strang
# put nucleotide strangs in list in one cell


input = []
for o in list1:
    if o.endswith('.fastq_AlignmentReport.txt'):
        Alignment = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+o, sep=':')
        avalues = (Alignment.iloc[5,1]) 
        input.append([o[0:o.index('.fastq_AlignmentReport.txt')],avalues])
#print (input)

df_input = pd.DataFrame(input)
df_input.columns = ['Experiment', 'Input']
df_output = pd.DataFrame(output)
df_output.columns = ['Experiment', 'Output']
#print (df_input)
#print (df_output)

# Merging the inputs and outputs
df_combined = pd.merge(df_input, df_output, on='Experiment', sort=True)
#print (df_combined)

isinstance(df_combined, pd.DataFrame) # is df_combined a pandas dateframe?  Yes. I checked because something wasnt't working as expected

#IMPORTANT - IF DATA IS NOT STORED AS float YOU CANNOT WORK WITH IT MATHEMATICALLY
df_combined['Input']=df_combined['Input'].astype(float)
df_combined['Output']=df_combined['Output'].astype(float)
# (df_combined.dtypes)

# Input columns substracte by Output columns and appended into a new columns
df_combined['seqslost'] = df_combined['Input'] - df_combined['Output'] # loss of sequences
#print (df_combined)

# Seqslost graph without stacking
# df_combined.plot.bar(x='Experiment', y=['Input', 'Output'], figsize=(10,3), ylabel='CloneCount')

# Seqslost graph with stacking
# df_combined.plot.bar(x='Experiment', y=['Output', 'seqslost'], stacked=True, figsize=(10,3), rot=90, ylabel='Total_Reads')

# With this script I have created a new dataframe (delete_ones) 
# In this dataframe, all the rows with cloneCounts of 1 are removed
# From here we have then extracted the new cloneCounts
# Ask Puneet about function, i am unsure how it works link below
# https://stackoverflow.com/questions/63280249/pandas-remove-rows-with-value-lower-than-a-threshold-but-keep-nans
output_del_1 = []
for b in list1:
    if b.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
        df_temp = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+b, sep='\t')
        df_temp = df_temp[(df_temp['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
        delete_ones = df_temp[(df_temp['cloneCount'] > 1)] # remove rows which have counts lower than 2
        del_1_values = (str(delete_ones.cloneCount.sum()))
        output_del_1.append([b[0:b.index('_UniqueCDR3_Exp_UniqueCDR3.txt')], del_1_values])
# print (output_del_1)


# Next we will merge these counts to the other dataset (df_combined) in a new dataframe called (df_combined2)
# Merge is used because we need to be able to connect them to their correct library
# Before merge, we will put them in list format and label the columns (there might be a smarter way but not for me)

df_output_del_1 = pd.DataFrame(output_del_1)
df_output_del_1.columns = ['Experiment', 'Output_del_ones']
df_output_del_1
df_combined = pd.merge(df_combined, df_output_del_1, on='Experiment', sort=True)
df_combined['Output_del_ones'] = df_combined['Output_del_ones'].astype(float)
# print (df_combined2)

df_combined['Samples_with_1_read'] = df_combined['Output'] - df_combined['Output_del_ones']
#print (df_combined2)



#df_combined2.plot.bar(x='Experiment', y=['Output', 'seqslost', 'Samples_with_1_read'], stacked=True, figsize=(10,3), rot=90, ylabel='CloneCounts')


cloneIds_count = []
for c in list1:
    if c.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
        df_temp = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+c, sep='\t')
        df_temp = df_temp[(df_temp['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
        temp1 = df_temp[(df_temp['cloneCount'] > 1)]
        temp2 = (str(df_temp.cloneId.count())) # count rows with the rows with 1 read included
        temp3 = (str(temp1.cloneId.count())) # count rows without the rows with 1 read included
        cloneIds_count.append([c[0:c.index('_UniqueCDR3_Exp_UniqueCDR3.txt')], temp2, temp3])
#print (cloneIds_count)

df_output_counted = pd.DataFrame(cloneIds_count)
df_output_counted.columns = ['Experiment', 'Uniqueclones_>0', 'Uniqueclones_>1']
df_combined = pd.merge(df_combined, df_output_counted, on='Experiment', sort=True)
#print (df_combined3)

#convert values in the following collumns to float values
df_combined['Uniqueclones_>0'] = df_combined['Uniqueclones_>0'].astype(float)
df_combined['Uniqueclones_>1'] = df_combined['Uniqueclones_>1'].astype(float)
#introduce a two digit labelling i.e. F1 = F01, and sort df_combined3 again
# dont include in pipeline
df_combined = df_combined.replace(regex=['F1_2', 'F2_2', 'F3_2', 'F5_2', 'F6_2', 'F7_2', 'F8_2', 'F9_2'], value=['F01_2', 'F02_2', 'F03_2', 'F05_2', 'F06_2', 'F07_2', 'F08_2', 'F09_2'])
df_combined = df_combined.sort_values(by='Experiment')
#print (df_combined3)
df_combined.to_csv('df_combined3.csv')
# ax = df_combined3.plot.bar(x='Experiment', y=['Uniqueclones_>1', 'Samples_with_1_read'], stacked=True, figsize=(10,3), rot=90, ylabel='Unique Clones')
# ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()]) #add thousand seperator

#%% Over or undersampling plots

#sequenced_deep_enough.py should be here

#%% Matching values between datasets - 2 hour run time - nucleotides
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
                Export = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
                Export = Export[(Export['cloneCount'] > 1)]   # by running this i get rid of the clones with only 1 cloneCount
                nt_sequences = (Export.nSeqCDR3)
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
        
#finaldataframe.to_excel("finaldataframe.xlsx") 

#%% converting F1 -> F01 etc - not sure which way is best
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
#%% Sorting in projects potentially - Should be changed to new datasets
#back_up_df = pd.read_excel('/Users/chvis/Jupyter-Demo/pandas_test/finaldataframe.xlsx')
#Now i need to use loc or iloc to slice out the data i need to use - each project alone
# i removed
Chris_data = back_up_df.iloc[[1,2,3,4,5,6,7,8,9,10,33,34], [1, 3,4,5,6,7,8,9,10,11,12,35,36]]
Line_data = back_up_df.iloc[[11,12,13,14,15,16,17,18,19,20,21,0], [1,13,14,15,16,17,18,19,20,21,22,23,2]]
Helen_data = back_up_df.iloc[[22,23,24,25,26,27,28,29], [1, 24,25,26,27,28,29,30,31]]
Myo2_data = back_up_df.iloc[[30,31,32], [1,32,33,34]]
Isabel_data = back_up_df.iloc[[35,36,37], [1, 37,38,39]]
Shirin_data = back_up_df.iloc [[38,39,40],[1,40,41,42]]
back_up_df_sortedinprojects = back_up_df.iloc[[1,2,3,4,5,6,7,8,9,10,33,34,11,12,13,14,15,16,17,18,19,20,21,0,22,23,24,25,26,27,28,29,30,31,32,35,36,37,38,39,40],[1,3,4,5,6,7,8,9,10,11,12,35,36,13,14,15,16,17,18,19,20,21,22,23,2,24,25,26,27,28,29,30,31,32,33,34,37,38,39,40,41,42]]
#Chris_data = back_up_df.loc[[0,[2:4],7], :]
#Chris_data
back_up_df_sortedinprojects

#%% Matching values between datsets for only top 100 clones - nucleotides - Nucleotide_match_tags.csv
## Copy of above, but mofified to only be top 100 clones
## This script takes 2 hours to run, therefore it has been commented out
## When you want to run it, remove the # from the four lines below here
## This script extract the nucleotide sequences and counts how many matches are between them and different samples
finaldataframe_top100 = pd.DataFrame()
addedlibcolumn = False #Checks whether the "Libs" column has been filled
libsnamestoadd = [] #Used to add the libraries to the "Libs" column
for Match_from_list in list1:
    arraydata = []
    dataframedata = []
    if addedlibcolumn == False:
        finaldataframe_top100['Libs'] = libsnamestoadd

    if Match_from_list.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
        Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+Match_from_list, sep='\t')
        Export = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
        Export = Export.iloc[0:100,0:]   # by running this i get rid of the clones with only 1 cloneCount
        nt_sequences = (Export.nSeqCDR3)
        df_Standard = pd.DataFrame(nt_sequences) #Make into actual dataframe - can be checked with ifinstance funciton
        
        for Match_to_list in list1:
                
            if Match_to_list.endswith('_UniqueCDR3_Exp_UniqueCDR3.txt'):
                if addedlibcolumn == False:
                     libsnamestoadd.append(Match_to_list[0:Match_to_list.index('_UniqueCDR3_Exp_UniqueCDR3.txt')])
                print([Match_to_list[0:Match_to_list.index('_UniqueCDR3_Exp_UniqueCDR3.txt')]])
                Export = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/'+Match_to_list, sep='\t')
                Export = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)               
                Export = Export.iloc[0:100,0:]   # by running this i get rid of the clones with only 1 cloneCount
                nt_sequences = (Export.nSeqCDR3)
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
            finaldataframe_top100['Libs'] = libsnamestoadd # Filles the empty 'Libs'column with the library names
            addedlibcolumn = True
        
        tempdataframe = pd.DataFrame(dataframedata)
        finaldataframe_top100[Match_from_list[0:Match_from_list.index('_UniqueCDR3_Exp_UniqueCDR3.txt')]] = tempdataframe
        print(finaldataframe_top100)
        #finaldataframe = pd.merge(finaldataframe, tempdataframe, on=Match_from_list)
finaldataframe_top100.to_csv('finaldataframe_top100.csv')
#%% converting and sorting full nucleotide dataset - post remove 1's and nt seqs not divisible by 3
os.chdir('/Users/chvis/Jupyter-Demo/pandas_test') 
Nucleotide_match = pd.read_excel('finaldataframe.xlsx',index_col=(1)).drop(labels="Unnamed: 0", axis=1)
convert_F1_F01_etc(Nucleotide_match)
Nucleotide_match = Nucleotide_match.sort_index(axis=0)
Nucleotide_match = Nucleotide_match.sort_index(axis=1)
tags = pd.read_excel('/Users/chvis/Jupyter-Demo/pandas_test/Experiments_extra_information.xlsx',index_col=(0))

Nucleotide_match_tags = tags.join(Nucleotide_match)
Nucleotide_match_tags.to_csv('Nucleotide_match_tags.csv')
#%% converting and sorting top100
finaldataframe_top100_sorted = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/finaldataframe_top100.csv',index_col=(1)).drop(labels="Unnamed: 0", axis=1)
convert_F1_F01_etc(finaldataframe_top100_sorted)
finaldataframe_top100_sorted = finaldataframe_top100_sorted.sort_index(axis=0)
finaldataframe_top100_sorted = finaldataframe_top100_sorted.sort_index(axis=1)
#finaldataframe_top100_sorted.to_csv('finaldataframe_top100_sorted.csv')   # only uncomment this when you want to export
print (finaldataframe_top100_sorted)

tags = pd.read_excel('/Users/chvis/Jupyter-Demo/pandas_test/Experiments_extra_information.xlsx',index_col=(0))

#join the matrix with the file containing library labels and project labels 
top100_sort_tags = tags.join(finaldataframe_top100_sorted)
top100_sort_tags.to_csv('top100_sort_tags.csv')
finaldataframe_top100_sorted.to_csv('finaldataframe_top100_sorted.csv')
#%% heat map get puneet to explain the "for i inxxx and for j in xxx"
#Need to figure out how to do the library and project tags in python instead of excel

os.chdir('/Users/chvis/Jupyter-Demo/pandas_test') 
df_heat = df_combined3[:]
df2_heat = finaldataframe_top100_sorted[:]
df3_heat = np.zeros([41,41]) # is the final matrix
for i in range(len(df3_heat)): # loop for having x and y axis
    for j in range(len(df3_heat)):
        df3_heat[i,j]= df2_heat.iloc[i,j]/((df_heat.loc[df_heat['Experiment']==df2_heat.index[i],
                                                        "Uniqueclones_>1"].iloc[0]+df_heat.loc[df_heat['Experiment']==df2_heat.columns[j],
                                                        "Uniqueclones_>1"].iloc[0])/2)
        # df.query('Experiment=='+df2.index[j])["Uniqueclones_>1"]
        # df3[i,j]= 1
#df2.columns
# df["Uniqueclones_>1"]
df4_heat = pd.DataFrame(df3_heat, columns= df2_heat.columns, index= df2_heat.index)
#df4_heat.to_csv('df4_heat.csv')

# in the 6 lines below i try to make a new column that has values of 1-4 based on library name yo use for library color labels
# However, i ended up giving up and doing it in excel after a few hours because i was in a rush - problem for future chris
#df5 = df4[:]
#df5.insert(0 , 'new_col', 'nan')
#for x in df5.index:
    #for y in df5.new_col:
        #if x.startswith('Library_1'):
            #df5.replace(to_replace=y, value=3, inplace=True)
#df5        

#plt.imshow(df4)
#plt.savefig("heatmap.png")
# df3.columns = df2.columns
# df4.index = df3.index

sns.heatmap(df4_heat, cmap="Blues", annot=True, annot_kws={"size":4}, fmt='.2f')
fig, ax = plt.subplots(figsize=(15,12), dpi = 400) 

#Color the samples based on libraries and projects
df_cluster = pd.read_csv('df4_heat_Lib_tags.csv', sep=';', index_col=(0))

Lib_color = df_cluster.Lib_Tag_y.map({
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

numerical_only = df_cluster.columns[2:]
modify = sns.clustermap(df_cluster[numerical_only], figsize=(12,12), annot=True, annot_kws={"size":4}, fmt='.2f', row_colors=[Lib_color, Project_color])
#plt.savefig('clustermap_L_P.png', dpi=400)
#modify.ax_row_dendrogram.set_visible(False)
### HERE IS THE ISSUE
#modify.ax_heatmap.xaxis.set_label_position('top')
#modify.ax_heatmap.yaxis.set_label_position('left')
###
#%% Don't use, just in case something goes wrong - Puneets unmodified script for heatmap

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created on Fri Jan  7 15:54:46 2022

#@author: puneetr
#"""
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

df = pd.read_csv("df_combined3_.csv")
df2 = pd.read_csv("back_up_df_sortedinprojects_.csv",index_col=(1)).drop(labels="Unnamed: 0", axis=1)
df3 = np.zeros([41,41])
for i in range(len(df3)):
    for j in range(len(df3)):
        df3[i,j]= df2.iloc[i,j]/((df.loc[df['Experiment']==df2.index[i], "Uniqueclones_>1"].iloc[0]+df.loc[df['Experiment']==df2.columns[j], "Uniqueclones_>1"].iloc[0])/2)
        # df.query('Experiment=='+df2.index[j])["Uniqueclones_>1"]
        # df3[i,j]= 1
#df2.columns
# df["Uniqueclones_>1"]
#df4 = pd.DataFrame(df3, columns= df2.columns, index= df2.index)
#plt.imshow(df4)
#plt.savefig("heatmap.png")
# df3.columns = df2.columns
# df4.index = df3.index


