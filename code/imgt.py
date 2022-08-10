#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 17:23:57 2022

@author: chvis
"""

# imgt code - not done

imgt = pd.read_csv('/Users/chvis/Jupyter-Demo/pandas_test/2_IMGT-gapped-nt-sequences.txt', sep='\t')
#imgt.columns
CDR3s_IMGT = (imgt.loc[:,'CDR3-IMGT'])
#CDR3s_IMGT
CDR3s_IMGT

#avalues = (Alignment.iloc[5,1]) 

imgt_count = imgt.pivot_table(columns=['CDR3-IMGT'], aggfunc='size')
print (imgt_count)
imgt_count.sort_values(ascending=False)