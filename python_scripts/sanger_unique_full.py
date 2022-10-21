#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 17:00:42 2022

@author: chvis
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np

## For each sanger sample, get unique sequences and extract their CDR3s - this means that the CDR3H can be identical


## Use this to get the full sequences combined
os.chdir('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/')
Sanger_combined = pd.read_csv('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/5_AA-sequences_scFv.txt', sep='\t').drop(labels="Unnamed: 18", axis=1)
## heavy chain
Full_H_info = Sanger_combined.dropna(subset=['V-D-J-REGION']) ## H dataset is equal to full dataset but the L-rows do not have the VDJ REGION and therefore this will remove them from this dataset
Full_H_info = Full_H_info.drop(['V-J-REGION', 'V-DOMAIN analysis order'], axis=1) ## drop the now empty columns
Full_H_info['V-DOMAIN ID'] = Full_H_info['V-DOMAIN ID'].str[:-2] ## remove the _H,_K,_L to easier join them
## light chain
Full_L_info = Sanger_combined.dropna(subset=['V-J-REGION']) #### H dataset is equal to full dataset but the H-rows do not have the VJ REGION and therefore this will remove them from this dataset
Full_L_info = Full_L_info.drop(['V-D-J-REGION', 'D-GENE and allele', 'V-DOMAIN analysis order'], axis=1)  ## drop the now empty columns
Full_L_info['V-DOMAIN ID'] = Full_L_info['V-DOMAIN ID'].str[:-2] ## remove the _H,_K,_L to easier join them

## combine heavy and light
Full_HL_combined = Full_H_info.merge(Full_L_info, on='V-DOMAIN ID', suffixes=('_H', '_L')) 
Full_HL_combined = Full_HL_combined.rename(columns={'V-DOMAIN ID': 'TPL_ID'})
Full_HL_combined['TPL_ID'] = Full_HL_combined['TPL_ID'].str.replace ('_TPL','TPL') ## deletes the _ that are in front of some seqs (macrogen i believe)
Full_HL_combined['TPL_ID'] = [x[:14] for x in Full_HL_combined['TPL_ID']] ## (only take the 14 first characters of each ID, thereby we remove the primer info)
## this merge results in less than either dataset combined, so i assume that this is when the light or heavy chain is only present, which is a failed sanger sequence anyways

## written by Punnet - this corrects the A1 top A01 in the TPL0039 samples. If the last character is _ then insert 0 inbetween A and 1 etc -> A01
Full_HL_combined['TPL_ID'] = [str(x[:-2]+'0'+x[-2]) if x[-1]=='_' else x for x in Full_HL_combined['TPL_ID']]

## combined into full scFv sequences
Full_HL_combined['Full_AA_Seq'] = Full_HL_combined['V-D-J-REGION']+'GGGGSGGGGSGGGAS'+Full_HL_combined['V-J-REGION']

sanger_extract_all = Full_HL_combined.loc[:,['TPL_ID','JUNCTION_H','Full_AA_Seq']]
sanger_extract_all.to_csv('sanger_all_sequences.csv')
sanger_extract_unique = sanger_extract_all.drop_duplicates('Full_AA_Seq')
sanger_extract_unique = sanger_extract_unique[(sanger_extract_unique.Full_AA_Seq.str.contains(r'[X]') == False)] ## remove seqs with X in it = sequencing error
sanger_extract_unique.to_csv('sanger_all_unique_sequences.csv')


## get unique sequences on full sequences, split into Vl and VH
sanger_extract_Biogenity = Full_HL_combined.loc[:,['TPL_ID','V-D-J-REGION','V-J-REGION','Full_AA_Seq']]
sanger_extract_Biogenity = sanger_extract_Biogenity.drop_duplicates('Full_AA_Seq')
sanger_extract_Biogenity = sanger_extract_Biogenity[(sanger_extract_Biogenity.Full_AA_Seq.str.contains(r'[X]') == False)] ## remove seqs with X in it = sequencing error
sanger_extract_Biogenity = sanger_extract_Biogenity[(sanger_extract_Biogenity.Full_AA_Seq.str.contains(r'[*]') == False)] ## remove seqs with X in it = sequencing error
sanger_extract_Biogenity = sanger_extract_Biogenity.loc[:,['V-D-J-REGION','V-J-REGION']]
sanger_extract_Biogenity.to_csv('sanger_extract_Biogenity.csv')

duplicate_count = sanger_extract_unique.pivot_table(columns=['JUNCTION_H'], aggfunc='size') # count duplicates
duplicate_count_df = pd.DataFrame(duplicate_count)
sanger_extract_CDR3 = sanger_extract_unique.drop_duplicates('JUNCTION_H') # drop CDR duplicates
sanger_extract_CDR3 = sanger_extract_CDR3.merge(duplicate_count_df, on='JUNCTION_H')
sanger_extract_CDR3.rename(columns={0:'Duplicates'}, errors="raise", inplace=True)
sanger_extract_CDR3.to_csv('sanger_all_unique_CDR3s.csv')
