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
    essential = None
    nonessential = None
    test = None 
    
    # other parameters, optional to change 
    metric = 'b'
    threshold = 0.5
    weight_e = 2
    weight_n = 0
    
    # define later 
    p = None
    score = None 
    essential_misses = None
    essential_hits = None 
    
    def __init__(self, test, essential, nonessential, 
                 metric = 'b',threshold=0.5,weight_e=2, weight_n=0.0):
        '''
        Initializes the class instance, takes pandas dataframes. 
        test = dataframe with the list of log2fold of TEST vs. PRIMARY for each gene 
        essential = dataframe of essential genes determined via pipeline 
        nonessential = dataframe of nonessential genes determined via pipeline 
        metric = 'b' or 'g' for binary or graded matrix for calculating hit and miss scores  
            default = 'b' 
        threshold = threshold for log2fold change to be considered differentially expressed
            default = 0.5 
        weight_e = weight for essential gene
            default = 2
        weight_n = weight for nonessential gene 
            default = 0.5
        '''
        self.test = test 
        self.essential = essential
        self.nonessential = nonessential
        self.metric = metric 
        self.threshold = threshold
        self.weight_e = weight_e
        self.weight_n = weight_n
    
        
    # gets score based on metric initialized during object initialization  
    def score(self):
        if(self.metric == 'b'):
            return {'Score': self.__binary(self.test), 'p':self.__get_p()}
        
        else:
            return {'Score': self.__graded(self.test), 'p':self.__get_p()}
        
    # local function for binary weighting 
    def __binary(self, test):
        ''' 
        Returns the differentiation score for the initiated class. 
        Defines alignment based on global threshold for log2foldchange and a p<0.05 
        Treats each "hit/miss" as the same no matter the log2 fold change/FDR 
        '''
        
        diff = test[(abs(test['log2FoldChange']) >= self.threshold) & 
                            (test['pvalue'] <= 0.05)]
        
        aligned = test[(abs(test['log2FoldChange']) < self.threshold) |
                            (test['pvalue'] > 0.05)]
        
        # get hits and misses of essential and nonessential genes 
        hits_e = aligned[aligned['genes'].isin(self.essential['genes'])]
        hits_n = aligned[aligned['genes'].isin(self.nonessential['genes'])]

        miss_e = diff[diff['genes'].isin(self.essential['genes'])]
        miss_n = diff[diff['genes'].isin(self.nonessential['genes'])]
        
        
        hits = len(hits_e)*self.weight_e + len(hits_n)*self.weight_n
        # do we want to weight same as hits, or over-weight if there is a miss of an essential gene?
        misses = len(miss_e)*self.weight_e - len(miss_n)*self.weight_n
        
        # store as global variables 
        self.essential_hits = hits_e
        self.essential_misses = miss_e 
        
        # store essential gene hits/misses as global variables 
        return hits / misses 
        
   

    # local function for graded weighting 
    def __graded(self,test):
        ''' 
        Returns the differentiation score for the initiated class. 
        Defines alignment based on global threshold for log2foldchange and a p<0.05 
        Treats each "hit/miss" differently depending on the log2fold change and FDR of the gene 
        '''
        diff = test[(abs(test['log2FoldChange']) >= self.threshold) & 
                           (test['pvalue'] <= 0.05)]
       
        aligned = test[(abs(test['log2FoldChange']) < self.threshold) |
                           (test['pvalue'] > 0.05)]
        
        # get hits and misses of essential and nonessential genes 
        hits_e = aligned[aligned['genes'].isin(self.essential['genes'])]
        hits_n = aligned[aligned['genes'].isin(self.nonessential['genes'])]

        miss_e = diff[diff['genes'].isin(self.essential['genes'])]
        miss_n = diff[diff['genes'].isin(self.nonessential['genes'])]
        
        
        hits = 0
        misses = 0
        for gene in hits_e:
            hits += gene['log2FoldChange']*(1 + np.log10(gene['pvalue']))*self.weight_e
        for gene in hits_n:
            hits += gene['log2FoldChange']*(1 + np.log10(gene['pvalue']))*self.weight_n
        for gene in miss_e:
            misses += gene['log2FoldChange']*(np.log10(gene['pvalue']))*self.weight_e
        for gene in miss_n:
            misses += gene['log2FoldChange']*(np.log10(gene['pvalue']))*self.weight_n
        
        
        # store essential gene hits/misses as global variables 
        self.essential_hits = hits_e
        self.essential_misses = miss_e 
        
        return hits / misses 
    
    
    
    def __get_p(self):
        '''
        Runs 1000 trials, permutes gene names across whole dataset. 
        '''
        vals = np.zeros(1000)
        for i in len(0,1000):
            vals[i] = self.__run_null()
    
        # only get values higher than score 
        more_extreme = [x for x in vals if x >= self.score] 
        toRet = len(more_extreme)/len(vals)
        
        if(toRet == 0 or toRet == 1): 
            self.p = 0.001
        else:
            self.p = toRet 
    
    
    def __run_null(self):
        '''
        Shuffles the gene names of the entire dataset, then runs the scoring algorithm on the shuffled data.
        Returns the differentiation score of the randomly shuffled dataframe. 
        '''
        
        copy = self.test.copy()
        genes = sorted(self.test['genes'], key=lambda k: random.random()) # shuffle 
        copy['genes'] = genes
        if(self.metric == 'b'):
            return self.__binary(copy)
        else:
            return self.__graded(copy)
        
        
if __name__ == '__main__':
    # must run get_essentail_genes.py first 
    essential = pd.read_table('../data/essential_CM.tsv')
    nonessential = pd.read_table('../data/unessential_CM.tsv')
    
    
    

    
    
               
    
 
    
    
    
    
    
        
    




    
    
    