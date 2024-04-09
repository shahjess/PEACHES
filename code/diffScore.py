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
        self.genes_to_test = self.test[self.test['Gene'].isin(self.gene_list['Gene'])]
   
        
    # gets score based on metric initialized during object initialization  
    def score(self):
        self.score = self.__score(self.test)
        return {'Score': self.score, 'p':self.__get_p()}
        
    
        
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
            gene_test = test[test['Gene']==gene]
            # in case gene is not in test case 
            #print(gene_test)
            if(len(gene_test)>0):
                total += 1/max(1,abs(gene_test['log2FoldChange'].to_numpy()[0]))
                n+=1
            
        return total/n
        
        
    
    def __get_p(self):
        '''
        Runs 1000 trials, permutes gene names across whole dataset. 
        '''
        vals = np.zeros(100)
        
        for i in range(0,100):
            vals[i] = self.__run_null()
    
        # only get values higher than score 
        more_extreme = [x for x in vals if x >= self.score] 
        toRet = len(more_extreme)/len(vals)
    
        
        if(abs(toRet - 0.0) < 1e-9): 
            self.p = 0.01
        else:
            self.p = toRet 
        
        return self.p
    
    def __float_equals(float1, float2, tolerance=1e-9):
        return abs(float1 - float2) < tolerance
    
    def __run_null(self):
        '''
        Shuffles the gene names of the entire dataset, then runs the scoring algorithm on the shuffled data.
        Returns the differentiation score of the randomly shuffled dataframe. 
        '''
        
        copy = self.test.copy()
        genes = sorted(self.test['Gene'], key=lambda k: random.random()) # shuffle 
        copy['Gene'] = genes
        
        return self.__score(copy)
        
        
        
if __name__ == '__main__':
    # must run get_essential_genes.py first 
    #test = pd.read_table('../data/log_output_fibroblast.tsv',sep='\t')
    test = pd.read_table('../data/log_output_cm.tsv',sep='\t')
    #test = pd.read_table('../data/log_output_lung.tsv',sep='\t')
    # needs to be the same format! 
    test.rename(columns={'Gene Name': 'Gene'}, inplace=True)
    
   #  test = pd.read_table('../data/primary_goodDiff.tsv',sep='\t')
    gene_list = pd.read_table('../data/essential_CM.tsv',sep='\t')
    gene_list = gene_list.iloc[:,1:]
    
    
    s = diffScore(test,gene_list)
    
    score = s.score()
    print(score)
    

    
    
    
    

    
    
               
    
 
    
    
    
    
    
        
    




    
    
    