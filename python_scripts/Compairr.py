#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 17:24:38 2022

@author: chvis
"""

### Now to get the data ready for compairr
## installed using brew install torognes/bioinf/compairr

## These scripts are mainly to convert the data into a format that fits the compairr command-line tool


import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns

## to test i have just exported two files - have not tried exporting loop yet

os.chdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/')
files = [x for x in os.listdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/') if x.endswith("Library_1_F1_2_UniqueCDR3_Exp_UniqueCDR3.txt")] # will read all txt files. so do not keep any not required txt files
temp = pd.DataFrame()
for i in files:
    temp = pd.read_csv(i,sep = "\t")
    temp = temp[(temp['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
    temp = temp[(temp['cloneCount'] > 1)]   # by running this i get rid of the clones with only x cloneCount
    temp = (temp[(temp.aaSeqCDR3.str.contains(r'[*_]') == False)]) ## Removes the rows with AA sequence containing asterisk
    temp = temp.head(1000)
    temp = temp.rename(columns={'cloneId':'sequence_id', 'cloneCount':'duplicate_count', 'bestVGene':'v_call', 'bestJGene':'j_call', 'nSeqCDR3':'junction', 'aaSeqCDR3':'junction_aa',})
    temp = temp.drop(['cloneFraction', 'lengthOfCDR3', 'meanQualCDR3'], axis=1)
    column_names= ['sequence_id', 'duplicate_count', 'v_call', 'j_call', 'junction', 'junction_aa',]
    temp = temp.reindex(columns = column_names)
    temp['duplicate_count'] = temp['duplicate_count'].astype(int)
    temp.insert(0, 'repertoire_id', i[0:i.index('_2_UniqueCDR3_Exp_UniqueCDR3.txt')] )
    i = i[0:i.index('_2_UniqueCDR3_Exp_UniqueCDR3.txt')]
    temp.to_csv(i + '.tsv', sep="\t")
    


## For using compairr like Victor have proposed, i need to run the -x for each library against the all list - then the length of each repetoire list...
## ... will be the matches that list had with its counterpart


    
### Pull all sequences into ones dataset

os.chdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/')
files = [x for x in os.listdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/') if x.endswith("_UniqueCDR3_Exp_UniqueCDR3.txt")] # will read all txt files. so do not keep any not required txt files
temp = pd.DataFrame()
alltemp = pd.DataFrame()
for i in files:
    temp = pd.read_csv(i,sep = "\t")
    temp = temp[(temp['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
    temp = temp[(temp['cloneCount'] > 1)]   # by running this i get rid of the clones with only x cloneCount
    temp = (temp[(temp.aaSeqCDR3.str.contains(r'[*_]') == False)]) ## Removes the rows with AA sequence containing asterisk
    temp = temp.head(1000)  
    temp = temp.rename(columns={'cloneId':'sequence_id', 'cloneCount':'duplicate_count', 'bestVGene':'v_call', 'bestJGene':'j_call', 'nSeqCDR3':'junction', 'aaSeqCDR3':'junction_aa',})
    temp = temp.drop(['cloneFraction', 'lengthOfCDR3', 'meanQualCDR3'], axis=1)
    column_names= ['sequence_id', 'duplicate_count', 'v_call', 'j_call', 'junction', 'junction_aa',]
    temp = temp.reindex(columns = column_names)
    temp['duplicate_count'] = temp['duplicate_count'].astype(int)
    temp.insert(0, 'repertoire_id', i[0:i.index('_2_UniqueCDR3_Exp_UniqueCDR3.txt')] )
    alltemp = alltemp.append(temp)
alltemp.to_csv('compairr_all_list_1000' + '.tsv', sep="\t")




##After having done the matrixes i want to make a heatmap from them

os.chdir('/Users/chvis/Jupyter-Demo/pandas_test/mixcr_analysis/')

matrix_d0 = pd.read_csv('matrix_1000_d0.tsv', sep='\t')



## if i want to only pull out some of the sequences
## i want to construct the violin plots both with matches and with the multiply thing - ...
## ... this way i will get the non-quantified and the quantified values

## run this in the terminal, changed the number after -d to changed the level of substitutions allowed
## compairr -m compairr_all_list.tsv compairr_all_list.tsv -u -d 0 -o matrix_1000_d0.tsv      
## 





















