#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 11:33:39 2022

@author: chvis
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


## get the columns you want from the combined_chris dataset
os.chdir('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/')
Chris_main = pd.read_csv('Combined_Chris.csv').drop(columns='Unnamed: 0')
Chris_main = Chris_main.loc[:,['TPL_ID','Full_AA_Seq','Selection_strategy','Ecarpholin_S_DDEL_5ug/mL','Myotoxin_II_DDEL_5ug/mL','PLA2_n.naja(Uniprot:P15445)_DDEL_5ug/mL','Ecarpholin_S_CDEL_100nM', 'Myotoxin_II_CDEL_100nM','PLA2_n.naja(Uniprot:P15445)_CDEL_100nM']]

Chris_main['Counts'] = Chris_main.groupby('Full_AA_Seq')['Ecarpholin_S_DDEL_5ug/mL','Myotoxin_II_DDEL_5ug/mL','PLA2_n.naja(Uniprot:P15445)_DDEL_5ug/mL'].transform('sum')

df = Chris_main.groupby('Full_AA_Seq')[]

os.chdir('/Users/chvis/Dropbox/PhD Project/Projects/Main project - How cross-reactive can an antibody be made/Sanger Sequencing/For_analysis_Oslo/Sanger_analysis/')
Line_project = pd.read_excel('Combined_Line.xlsx').drop(columns='Unnamed: 0')
















				