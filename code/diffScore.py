#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:35:22 2024

@author: tessfallon

Contains a class diffScore() which defines a quantitative measurement of how well 
a cell is differentiated based on a list of essential genes specific to the cell type. 

After running diffScore.score(), global variables define the score, the p-value, and a list 
of essential genes aligned and missing within the tested RNAseq data. 

"""

import pandas as pd
import numpy as np
import random

class diffScore:
    
    # initialize in constructor 
  
    test = None 
    gene_list = None
    genes_to_test = None 

    
    # define later 
    p = None
    score = None 

    
    def __init__(self, test, gene_list):
        '''
        Initializes the class instance, takes pandas dataframes. 
        test = dataframe with the list of log2fold of PRIMARY vs. TEST for each gene 
        gene_list = dataframe with genes determined to be important for CM function
        '''
        self.test = test 
        self.gene_list= gene_list
        self.genes_to_test = self.test[self.test['Gene Name'].isin(self.gene_list['Gene Name'])]
   
        
    # gets score based on metric initialized during object initialization  
    def score(self):
        self.score = self.__score(self.test)
        return self.score
        
    
        
    # local function for binary weighting 
    def __score(self, test):
        ''' 
        Returns the differentiation score for the initiated class. 
        Defines alignment based on global threshold for log2foldchange and a p<0.05 
        Treats each "hit/miss" as the same no matter the log2 fold change/FDR 
        '''
        
        total = 0
        n = 0
        for i in range(0,len(self.gene_list)):
            gene = self.gene_list.iloc[i, 0]
            gene_test = test[test['Gene Name']==gene]
            # in case gene is not in test case 
            #print(gene_test)
            if(len(gene_test)>0):
                total += 1/max(1,abs(gene_test['log2FoldChange'].to_numpy()[0]))
                n+=1
            
        return total/n
        
        
    
        
        
        
if __name__ == '__main__':
    # must run get_essential_genes.py first 
    import matplotlib.pyplot as plt
    import statistics as stats 
    
    i1 = pd.read_table('../data/comparison_output/ipsc1/ipsc1_log_output.tsv',sep='\t') 
    i2 = pd.read_table('../data/comparison_output/ipsc2/log_output.tsv',sep='\t') 
    i3 = pd.read_table('../data/comparison_output/ipsc3/ipsc3_log_output.tsv',sep='\t') 
    
    p1 = pd.read_table('../data/comparison_output/primary1/primary1_log_output.tsv',sep='\t') 
    p2 = pd.read_table('../data/comparison_output/primary2/primary2_log_output.tsv',sep='\t') 
    p3 = pd.read_table('../data/comparison_output/primary3/primary3_log_output.tsv',sep='\t') 
    
    l1 = pd.read_table('../data/comparison_output/lung1/lung1_log_output.tsv',sep='\t') 
    l2 = pd.read_table('../data/comparison_output/lung2/lung2_log_output.tsv',sep='\t') 
    l3 = pd.read_table('../data/comparison_output/lung3/lung3_log_output.tsv',sep='\t') 
    
    gene_list = pd.read_table('../data/essential_CM.tsv',sep='\t')
    gene_list = gene_list.iloc[:,1:]
    
    gene_list.rename(columns={'Gene': 'Gene Name'}, inplace=True)
    
    s_i1 = diffScore(i1, gene_list)
    s_i2 = diffScore(i2, gene_list)
    s_i3 = diffScore(i3, gene_list)
    
    s_p1 = diffScore(p1, gene_list)
    s_p2 = diffScore(p2, gene_list)
    s_p3 = diffScore(p3, gene_list)
    
    s_l1 = diffScore(l1, gene_list)
    s_l2 = diffScore(l2, gene_list)
    s_l3 = diffScore(l3, gene_list)
    
    ipscs = np.array([s_i1.score(), s_i2.score(),s_i3.score()])
    ps = np.array([s_p1.score(), s_p2.score(),s_p3.score()])
    ls = np.array([s_l1.score(), s_l2.score(),s_l3.score()])
    
    
    # Create a figure and axis
    fig, ax = plt.subplots()
    
    # Plot the data points as dots
    ax.plot(np.ones_like(ls), ls, 'bo', label='Lung')
    ax.plot(2 * np.ones_like(ps), ps, 'go', label='Primary Cultured')
    ax.plot(3 * np.ones_like(ipscs), ipscs, 'ro', label='iPSC-Derived')
    
     # Add labels and title
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Differentiation Score')
    ax.set_title('Differentiation Scores by Cell Type')
    
    # Set x-axis ticks
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['Lung', 'Primary Culture', 'iPSC-Derived'])
    ax.set_ylim(0,1)
    
    # Show the plot
    plt.show()
    
    plt.savefig('../Figures/score_comparisons.png')
   
    



    
    
    
    
    ''' intial tests 
    test_fib = pd.read_table('../data/log_output_fibroblast.tsv',sep='\t')
    test_cm = pd.read_table('../data/log_output_cm.tsv',sep='\t')
    test_lung = pd.read_table('../data/log_output_lung.tsv',sep='\t')
    # needs to be the same format! 
    test_fib.rename(columns={'Gene Name': 'Gene'}, inplace=True)
    test_cm.rename(columns={'Gene Name': 'Gene'}, inplace=True)
    test_lung.rename(columns={'Gene Name': 'Gene'}, inplace=True)
    
    gene_list = pd.read_table('../data/essential_CM.tsv',sep='\t')
    gene_list = gene_list.iloc[:,1:]
    
    
    s_fib = diffScore(test_fib,gene_list)
    score_fib = s_fib.score()
    s_cm = diffScore(test_cm,gene_list)
    score_cm = s_cm.score()
    s_lung = diffScore(test_lung,gene_list)
    score_lung = s_lung.score()
    
    scores = [s_lung.score, s_fib.score, s_cm.score]
    labs = ['Lung', 'Direct Conversion', 'iPSC Derived']
    
    # Plotting
    plt.bar(labs, scores)
    
    # Adding labels
    plt.xlabel('Cell Type')
    plt.ylabel('Differentiation Scores')
    plt.title('Differentiation Scores by Cell Type')

    # Rotating x-axis labels for better readability
    plt.xticks(rotation=45)
    plt.ylim((0,1))
    
    plt.savefig('../Figures/score_barGraph.png')
    
    plt.show()
    ''' 

    

    
    
    
    

    
    
               
    
 
    
    
    
    
    
        
    




    
    
    