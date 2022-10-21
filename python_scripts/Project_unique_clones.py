#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 14:52:44 2022

@author: chvis
"""

import os
import pandas as pd
import random
import matplotlib.pyplot as plt

## In this script i want to show a bar graph for the unique sequences found in each project

os.chdir('/Users/chvis/Jupyter-Demo/pandas_test')

dataset = pd.read_csv('df_combined3.csv').drop(columns='Unnamed: 0')
dataset = dataset.replace({'Experiment':Library_2_to_panning})
dataset = dataset.sort_values(by='Experiment', ignore_index=True)

Chris_df = dataset.iloc[[0,1,5,9,2,3,4,6,7,8,10,11],:]
Line_df = dataset.iloc[[15,21,16,19,22,24,17,18,20,23,25,26],:]
Helen_df = dataset.iloc[[27,28,32,29,30,31,33,34],:]
Chris_myo_df = dataset.iloc[[12,13,14],:]
Isabel_df = dataset.iloc[[35,36,37],:]
Shirin_df = dataset.iloc[[38,39,40],:]

## Function to quickly make graphs
## PROBLEM: Savefig does not work in the function for some reason, so will do them manually instead
## Turn the last line on or off depending on whether you want to export

def make_barplot_unique_seqs(df, outfile):
    plt.figure()
    labels = df['Experiment']
    data = df['Uniqueclones_>1']
    plt.bar(labels,data)
    plt.ylabel('Unique clones')
    plt.xlabel('Selection rounds')
    plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
    plt.draw()
    plt.tight_layout()
    #vlines = [0.5, 3.5] ## Chris
    #vlines = [1.5, 5.5] ## Line
    vlines = [0.5, 2.5] ## Helen

    for vl in vlines:
        plt.axvline(x=vl, linestyle='--', color='black')
    plt.savefig(outfile+'_uni_clones.png', dpi=300)


make_barplot_unique_seqs(Chris_df, 'Chris')
make_barplot_unique_seqs(Line_df, 'Line')
make_barplot_unique_seqs(Helen_df, 'Helen')
make_barplot_unique_seqs(Chris_myo_df, 'Chris_myo')
make_barplot_unique_seqs(Isabel_df, 'Isabel')
make_barplot_unique_seqs(Shirin_df, 'Shirin')





















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

Library_to_TPL_1digit = {'Library_1_F1_' : 'TPL0049',
             'Library_1_F2_' : 'TPL0050',
             'Library_1_F3_' : 'TPL0065',
             'Library_1_F5_' : 'TPL0066',
             'Library_1_F6_' : 'TPL0067',
             'Library_1_F7_' : 'TPL0068',
             'Library_1_F8_' : 'TPL0095',
             'Library_1_F9_' : 'TPL0096',
             'Library_1_F10_' : 'TPL0097',
             'Library_1_F11_' : 'TPL0098',
             'Library_1_F12_' : 'TPL0099',
             'Library_2_F1_' : 'TPL0101',
             'Library_2_F2_' : 'TPL0103',
             'Library_2_F3_' : 'TPL0109',
             'Library_2_F5_' : 'TPL0110',
             'Library_2_F6_' : 'TPL0111',
             'Library_2_F7_' : 'TPL0123',
             'Library_2_F8_' : 'TPL0124',
             'Library_2_F9_' : 'TPL0125',
             'Library_2_F10_' : 'TPL0126',
             'Library_2_F11_' : 'TPL0127',
             'Library_2_F12_' : 'TPL0128',
             'Library_3_F1_' : 'TPL0214',
             'Library_3_F2_' : 'TPL0220',
             'Library_3_F3_' : 'TPL0221',
             'Library_3_F5_' : 'TPL0227',
             'Library_3_F6_' : 'TPL0266',
             'Library_3_F7_' : 'TPL0228',
             'Library_3_F8_' : 'TPL0229',
             'Library_3_F9_' : 'TPL0230',
             'Library_3_F10_' : 'TPL0037',
             'Library_3_F11_' : 'TPL0038',
             'Library_3_F12_' : 'TPL0039',
             'Library_4_F1_' : 'TPL0129',
             'Library_4_F2_' : 'TPL0130',
             'Library_4_F3_' : 'TPL0204',
             'Library_4_F5_' : 'TPL0239',
             'Library_4_F6_' : 'TPL0256',
             'Library_4_F7_' : 'TPL0285',
             'Library_4_F8_' : 'TPL0301',
             'Library_4_F9_' : 'TPL0302',   
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
             'TPL0130' : 'ACB_K',
             'TPL0204' : 'H--_L',
             'TPL0239' : 'HH-_L',
             'TPL0256' : 'HHH_L',
             'TPL0285' : 'J--_L_100uM',
             'TPL0301' : 'JJ-_L_1nM',
             'TPL0302' : 'JJJ_L_200pM',   
    }
