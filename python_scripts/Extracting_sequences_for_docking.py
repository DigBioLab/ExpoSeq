#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 09:31:02 2022

@author: chvis
"""

import os
import 

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

It was## get unique sequences on full sequences, split into Vl and VH
For_modelling = Full_HL_combined.loc[:,['TPL_ID','V-D-J-REGION','V-J-REGION','Full_AA_Seq', 'JUNCTION_H']]
For_modelling = For_modelling.drop_duplicates('Full_AA_Seq')
For_modelling = For_modelling[(For_modelling.JUNCTION_H.str.contains(r'[X]') == False)] ## remove seqs with X in CDR3H = sequencing error
For_modelling = For_modelling[(For_modelling.JUNCTION_H.str.contains(r'[*]') == False)] ## remove seqs with * in CDR3H = sequencing error
For_modelling.Full_AA_Seq = For_modelling.Full_AA_Seq.str.replace('X', 'A') ## replaces X with A
For_modelling.Full_AA_Seq = For_modelling.Full_AA_Seq.str.replace('*', 'A') ## replaces X with A

For_modelling = For_modelling.loc[:,['TPL_ID','V-D-J-REGION','V-J-REGION','Full_AA_Seq']]
For_modelling = For_modelling.rename(columns={'V-D-J-REGION':'Heavy_chain', 'V-J-REGION':'Light_chain'})
For_modelling.to_csv('Sanger_seqs_for_modelling.csv')