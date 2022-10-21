#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 14:30:30 2022

@author: chvis
"""

import os
import pandas as pd

## import biogenity output
os.chdir('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/Biogenity/ABPred-Seqbased-results')
biogenity_output = pd.read_csv('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/Biogenity/ABPred-Seqbased-results/output.csv')
biogenity_output['Full_AA_Seq'] = biogenity_output['VH']+'GGGGSGGGGSGGGAS'+biogenity_output['VL']
biogenity_output = biogenity_output.drop(columns=['VH','VL'])


## import the binding per project datasets to merge with
## import by switching working directory because it seemed easier
## best to keep them seperate because i have not found a good way of showing them in one dataset with all their different signals
os.chdir('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis')

Chris = pd.read_csv('Chris_main_df_full.csv')
Line = pd.read_csv('Line_project_df_full.csv')
Chris_myo = pd.read_csv('Chris_myo_project_df_full.csv')
Helen = pd.read_csv('Helen_project_df_full.csv')

## Use biogenity folder as working directory
os.chdir('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/Biogenity')

## merge based on full sequences
Chris = Chris.merge(biogenity_output, on='Full_AA_Seq')
Line = Line.merge(biogenity_output, on='Full_AA_Seq')
Chris_myo = Chris_myo.merge(biogenity_output, on='Full_AA_Seq')
Helen = Helen.merge(biogenity_output, on='Full_AA_Seq')
