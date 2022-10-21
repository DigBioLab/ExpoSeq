#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 00:06:40 2022

@author: chvis
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sn


## Quick check if all the sanger seqs are in the all unique seq list
        #matches[j,k] = len([x for x in list_all[j] if x in unique[k]])
        
## all_seq_combined = all the sequences found in the NGS
## CDR3_only = sanger matches somewhat trimmed

CDR3_only['nSeqCDR3'] = CDR3_only['nSeqCDR3'].str.upper()
all_seq_df = pd.DataFrame(all_seq_combined)
[x for x in CDR3_only['nSeqCDR3'] if x not in all_seq_combined]


nonmatched = pd.DataFrame()
for x in CDR3_match['nSeqCDR3']:
    if x not in all_seq_combined:
        nonmatched = nonmatched.append(CDR3_match.loc[CDR3_match['nSeqCDR3'] == x])
        
        print()


## take out all Sanger sequences + TPL_IDs that does not match their respected NGS round - all Sanger sequences should be located - otherwise NGS does not make sense
## NEEDED - RUN SANGER_PROJECT_HEATMAP first 3 sections

all_seqs_in_groups_series
CDR3_match_grouped
### use these two!!

new_df = all_seqs_in_groups_series.merge(CDR3_match_grouped, on='left', how='outer')

for j in range(len(list_all)):
    for k in range(len(unique)):
        matches[j,k] = len([x for x in list_all[j] if x in unique[k]])

for x in matches_df:
    for y in matches_df:
        print(x)
        print(y)

CDR3_match = CDR3_only[:]
CDR3_match = CDR3_match.drop_duplicates(subset = ['TPL_ID']) # remove duplicate TPL_IDs, some of Helens were sequenced twice
CDR3_match = CDR3_match.dropna(subset=['nSeqCDR3'])
CDR3_match = CDR3_match.sort_values('TPL_ID')
CDR3_match['nSeqCDR3'] = CDR3_match['nSeqCDR3'].str.upper() # matching is case sensitive so we need to conver to same case 
group = CDR3_match.groupby(by=['TPL_ID_short'])
#group.size()
CDR3_match_grouped = CDR3_match.groupby(by=['TPL_ID_short'])['nSeqCDR3'].apply(list) # group by TPL_ID_short and make list with all seqs in it
group = CDR3_match.groupby(by=['TPL_ID_short'])['nSeqCDR3']
total = group.size()




os.chdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/')
files = [x for x in os.listdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/') if x.endswith("_UniqueCDR3_Exp_UniqueCDR3.txt")] # will read all txt files. so do not keep any not required txt files
list_all = []
name_all = [] 
for i in files:
    temp = pd.read_csv(i,sep = "\t")
    temp = temp[(temp['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
    temp = temp[(temp['cloneCount'] > 1)]   # by running this i get rid of the clones with only x cloneCount
    #list_all.append(','.join(temp.nSeqCDR3.to_list())) ### this might have been the code that slowed it down
    name_all.append(i[0:i.index('2_UniqueCDR3_Exp_UniqueCDR3.txt')]) # remove unwanted ending of names
all_seqs_in_groups = pd.DataFrame(list_all)
all_seqs_in_groups['TPL_ID_short'] = name_all

all_seqs_in_groups.rename(columns=Library_to_TPL_1digit, inplace=True)
new_df = all_seqs_in_groups.merge(CDR3_match_grouped, on='TPL_ID_short', how='outer')
new_df['TPL_ID_short'].rename(Library_to_TPL_1digit, inplace=True)





CDR3_match_group_df = pd.DataFrame(CDR3_match_grouped)

## removing duplicates

## get only the unique CDR3s for each selection round
unique_only = []
for x in CDR3_match_grouped:
    unique_only.append(list(set(x))) 
CDR3_match_group_df['unique_only'] = unique_only



## grab NGS mixcr analysis



CDR3_match_group_df = CDR3_match_group_df.drop(columns='nSeqCDR3')
CDR3_match_group_df.join()


##
##
## i want to compare the Sanger sequencing to its respective NGS sequencing - then from this i want to pull out the non-matches
##
##


##all_seqs_in_groups ## make into series
## then all_seqs_in_groups.index = name_all - once it is a series that is


test_df = pd.DataFrame()
for y in test_df:
    if all_seqs_in_groups_series.index == CDR3_match_group_df['TPL_ID_short']:
        test_df = test_df.append([x for x in CDR3_match_group_df['unique_only'] if x in all_seqs_in_groups_series])
        print('yes')
    


for j in range(len(all_seqs_in_groups)):
    for k in range(len(unique)):
        if k :
            
        [x for x in list_all[j] if x in unique[k]])



for x in CDR3_match_group_df['unique_only']:
    if CDR3_match_group_df['TPL_ID'] == all_seqs_in_groups_df[]
    print (x)
    



nonmatched = nonmatched.append(CDR3_match.loc[CDR3_match['nSeqCDR3'] == x])




## get the matches between the NGS samples and the Sanger samples by using the unique Sanger samples per TPL_label as calculated further up
matches = np.zeros([41,27]) ## change it to sample size
for j in range(len(list_all)):
    for k in range(len(unique)):
        matches[j,k] = len([x for x in list_all[j] if x in unique[k]])
        matches_df = pd.DataFrame(matches, index=name_all, columns=CDR3_match_grouped.index)
        matches_df.rename(index=Library_to_TPL_1digit, inplace=True)

        