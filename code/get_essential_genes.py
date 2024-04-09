#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 13:30:09 2024

@author: tessfallon

Script that reads the files from the first pipeline step as .tsv files with 
column names [Gene log2FoldChange -log10FDR Product]

Finds the subpopulation of genes that are cardiac-specific and necessary for good differentiation 


"""

import pandas as pd

'''
# run once to split up data into essential and nonessential genes, save as .tsv files  

# read in datafiles 
og_primary = pd.read_table('../data/DE_og_primary.tsv')
goodDiff_primary = pd.read_table('../data/DE_allDiff_primary.tsv')
goodDiff_badDiff = pd.read_table('../data/DE_bad_good.tsv')

# genes needed for good cardiac differentiation 
diff_good = goodDiff_badDiff[((goodDiff_badDiff['log2FoldChange'] >= 0.5) & 
                         (goodDiff_badDiff['-log10FDR'] >= 1.31))]




# genes not needed--> cutoff 2 or 5 depending on title of file 
unessential = goodDiff_primary[((goodDiff_primary['log2FoldChange'] <= -2) & 
                         (goodDiff_primary['-log10FDR'] >= 1.31))]

unessential =  unessential[unessential['Gene'].isin(diff_good['Gene'])]


# cm essential genes 
essential =  diff_good[diff_good['Gene'].isin(unessential['Gene']) == False]


essential.to_csv('../data/essential_CM_2.tsv',sep='\t')
unessential.to_csv('../data/unessential_CM_2.tsv',sep='\t')


'''

# read in datafiles 
og_primary = pd.read_table('../data/DE_og_primary_new.tsv')
goodDiff_primary = pd.read_table('../data/primary_goodDiff.tsv')
goodDiff_badDiff = pd.read_table('../data/DE_bad_good.tsv')

# CM-specific genes are differentially expressed in primary and not og 
cm_specific = og_primary[((og_primary['log2FoldChange'] >= 0.5) & 
                         (og_primary['-log10FDR'] >= 1.31))]   # 1.3 is -log10(0.05)

# genes needed for cardiac differentiation 
diff_good = goodDiff_badDiff[((goodDiff_badDiff['log2FoldChange'] >= 0.5) & 
                         (goodDiff_badDiff['-log10FDR'] >= 1.31))]


# genes needed for cardiac differentiations that are CM-specific 
cm_needed = cm_specific[cm_specific['Gene'].isin(diff_good['Gene'])]


# genes not needed
unessential = goodDiff_primary[((goodDiff_primary['log2FoldChange'] <= -2) & 
                         (goodDiff_primary['-log10FDR'] >= 1.31))]

unessential =  unessential[unessential['Gene'].isin(cm_needed['Gene'])]


# cm essential genes 
essential =  cm_needed[cm_needed['Gene'].isin(unessential['Gene']) == False]



essential.to_csv('../data/essential_CM.tsv',sep='\t')
unessential.to_csv('../data/unessential_CM.tsv',sep='\t')







