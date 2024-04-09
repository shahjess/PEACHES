#pip install pydeseq2
#%pip install scanpy
#%pip install sanbomics
#%pip install bioinfokit

import PySimpleGUI as sg
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from functools import reduce
from sanbomics.plots import volcano
from bioinfokit import analys,visuz
from functools import reduce


import pandas as pd
import seaborn as sns
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import sys
sys.setrecursionlimit(10000)  # Set a higher recursion limit (adjust the value as needed)
import gseapy as gp
from gseapy.plot import gseaplot

def main():
    # Set the theme
    sg.theme('DarkTeal12')

    # Define the layout of the GUI
    layout = [
        [sg.Text('Directory Location for Stem Cell-Derived, Reprogrammed, or Converted Cells:', size=(60, 1))],
        [sg.InputText(key='-IPSC_CELL-', size=(50, 1))],
        [sg.Text('Directory Location for Primary Cells:', size=(60, 1))],
        [sg.InputText(key='-PRIMARY_CELL-', size=(50, 1))],
        [sg.Button('Submit'), sg.Button('Cancel')]
    ]

    # Create the window
    window = sg.Window('Please Upload Your Transcriptome Datasets', layout, finalize=True)

    # Event loop
    while True:
        event, values = window.read()

        if event == sg.WINDOW_CLOSED or event == 'Cancel':
            break
        elif event == 'Submit':
            ipsc_cell_dir = values['-IPSC_CELL-']
            primary_cell_dir = values['-PRIMARY_CELL-']
            if ipsc_cell_dir == '' or primary_cell_dir == '':
                sg.popup_error('Please upload directory for iPSC and Primary Cells!')
            else:
                sg.popup('Transcriptome sets uploaded successfully! Now running PEACHES.')
                log_output = differential_expression(ipsc_cell_dir, primary_cell_dir)
                print(log_output)
                window.close()
def create_meta_data(ipsc_cell_dir, primary_cell_dir):
    # Initialize empty arrays to store dataframes
    ipsc_dataframes = []
    primary_dataframes = []
    print("creating meta_data...")
    # Iterate over files in the iPSC cell directory
    for file in os.listdir(ipsc_cell_dir):
        
        if file.endswith(".tsv"): 
          ipsc_dataframes.append(pd.read_table(os.path.join(ipsc_cell_dir, file), sep='\t', index_col=False))

    # Iterate over files in the primary cell directory
    for file in os.listdir(primary_cell_dir):
        if file.endswith(".tsv"):
            primary_dataframes.append(pd.read_table(os.path.join(primary_cell_dir, file), sep='\t', index_col=False))
  
    # List to store the modified and renamed DataFrames
    modified_dfs = []

    primary = 0
    esc = 0
    ipsc = 0
    undiff = 0
    bad = 0
    
    # Process iPSC dataframes
    for i, sample in enumerate(ipsc_dataframes):
        subset = sample[['Gene Name', 'Expression Data']].copy()
        sample_name = sample.loc[0,'Characterization']
        
        if(sample_name=='Primary' or sample_name=='primary' or sample_name=='uiPSC'):
            sample_name = 'p'+str(primary)
            primary+=1
        elif(sample_name=='ESC'):
            sample_name = 'e'+str(esc)
            esc+=1
        #elif(sample_name=='Undiff_ESC' or sample_name=='Undiff_iPSC' or sample_name=='uiPSC'):
            #sample_name = 'u'+str(undiff)
            #undiff+=1
            # only for good diff vs. bad diff metadata
            continue
        elif(sample_name == 'biPSC'):
            sample_name = 'b'+str(bad)
            bad+=1
        elif(sample_name == 'iPSC' or sample_name=='fiPSC'):
            sample_name = 'i'+str(ipsc)
            ipsc+=1
        else:
            sample_name = 'i'+str(ipsc)
            ipsc+=1
        subset.sort_values('Gene Name', inplace=True)
        subset = subset.rename(columns={'Expression Data': sample_name})
        modified_dfs.append(subset)

    # Process primary cell dataframes
    for i, sample in enumerate(primary_dataframes):
        subset = sample[['Gene Name', 'Expression Data']].copy()
        sample_name = sample.loc[0,'Characterization']
        
        if(sample_name=='Primary' or sample_name=='primary'or sample_name=='uiPSC'):
            sample_name = 'p'+str(primary)
            primary+=1
        elif(sample_name=='ESC'):
            sample_name = 'e'+str(esc)
            esc+=1
        #elif(sample_name=='Undiff_ESC' or sample_name=='Undiff_iPSC' or sample_name=='uiPSC'):
            #sample_name = 'u'+str(undiff)
            #undiff+=1
            # only for good diff vs. bad diff metadata
            continue
        elif(sample_name == 'biPSC'):
            sample_name = 'b'+str(bad)
            bad+=1
        elif(sample_name == 'iPSC'):
            sample_name = 'i'+str(ipsc)
            ipsc+=1
        elif(sample_name == 'fiPSC'):
            sample_name = 'i'+str(ipsc)
            ipsc+=1
        else:
            sample_name = 'p'+str(primary)
            primary+=1
        subset.sort_values('Gene Name', inplace=True)
        subset = subset.rename(columns={'Expression Data': sample_name})
        modified_dfs.append(subset)
    print(modified_dfs)
    # Merge all modified dataframes
    merged = reduce(lambda left, right: pd.merge(left, right, on='Gene Name', how='inner'), modified_dfs)
    print("h3")
    merged = merged.drop_duplicates(subset='Gene Name', keep='first')
    names = merged['Gene Name']
    merged = merged.drop('Gene Name',axis=1)

    copy = merged.T

    copy = copy.T
    copy['Gene Name'] = names
    copy = copy.set_index('Gene Name')
    copy = copy.T
     
    meta_data = copy
    
    return meta_data

def assign_label(cell):

    if str(cell).startswith('p'):
        return 'Primary'
    elif str(cell).startswith('u'):
        return 'Undiff_esc'
    elif str(cell).startswith('i'):
        return 'hiPSC'
    elif str(cell).startswith('b'):
        return 'biPSC'
    elif str(cell).startswith('e'):
        return 'esc'


def differential_expression(ipsc_cell_dir, primary_cell_dir):
    meta_data = create_meta_data(ipsc_cell_dir, primary_cell_dir)
    
    #run DeSEQ
    print("running DeSeq...")
    merged = meta_data
    
    merged.set_index(merged.index, inplace=True)
    
    merged.sort_index(inplace=True)
    
    merged = merged.astype(int)
    
    labels = []

    for cell in merged.index:

        labels.append(assign_label(cell))

    metadata = pd.DataFrame(zip(merged.index, labels),
                        columns = ['Sample', 'Condition'])

    metadata = metadata.set_index('Sample')
    metadata = metadata.dropna()
    print(metadata)
    dds = DeseqDataSet(counts=merged.iloc[0:20],
            metadata=metadata,
            design_factors="Condition")
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast = ('Condition','Primary','hiPSC'))

    stat_res.summary()
    res = stat_res.results_df
    stat_res_summary = stat_res.results_df
    #drop nan rows
    stat_res_summary = stat_res_summary.dropna()
    log_output = stat_res_summary
    
    desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    file_path = os.path.join(desktop_path, 'log_output.tsv')

    # Writing data to the TSV file
    log_output.to_csv(file_path, sep='\t')

    print("TSV log_fold file saved successfully on desktop.")

    #volcano
    res = log_output
    #Filter based on p-value
    ranking = res[res['pvalue'] <= 0.05]

    #Sort based on abs log2 fold change
    ranking = ranking[['log2FoldChange']].dropna()
    ranking['abs_l2FC'] = ranking['log2FoldChange'].abs()
    ranking.sort_values(by='abs_l2FC', ascending=False, inplace=True)
    ranking = ranking[ranking['abs_l2FC'] > 1]

    ranking.drop(columns=['abs_l2FC'], inplace=True)
    pre_res = gp.prerank(rnk=ranking, # or rnk = rnk,
                     gene_sets='GO_Biological_Process_2021',
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )
    out = []

    for term in list(pre_res.results):
        out.append([term,
                   pre_res.results[term]['fdr'],
                   pre_res.results[term]['es'],
                   pre_res.results[term]['nes'],
                   pre_res.results[term]['gene %'],
                   pre_res.results[term]['lead_genes'],])

    out_df = pd.DataFrame(out, columns = ['Term','fdr', 'es', 'nes','gene %','lead_genes']).sort_values('nes').reset_index(drop = True)

    out_df['abs_nes'] = abs(out_df['nes'])  # Create a new column with absolute NES values
    out_df = out_df.sort_values('nes', ascending=False).reset_index(drop=True)

    out_df.drop('abs_nes', axis=1, inplace=True)
    out_df = out_df.sort_values(by='nes', key=lambda x: abs(x), ascending=False)
    # Save the sorted DataFrame to a TSV file
    file_path_out = os.path.join(desktop_path, 'go_output.tsv')
    out_df.to_csv(file_path_out, sep='\t', index=False)
    
    out_df = out_df.reindex(out_df['nes'].abs().sort_values(ascending=False).index)
    primary_v_ipsc_sorted_top20 = out_df.head(20).copy()
    terms = primary_v_ipsc_sorted_top20['Term']
    nes = primary_v_ipsc_sorted_top20['nes'].abs() # Taking absolute value
  

    plt.figure(figsize=(10, 6))
    plt.barh(terms, nes, color='skyblue')
    plt.xlabel('Absolute Normalized Enrichment Score (|NES|)')  # x-axis label
    plt.ylabel('Gene GO Pathway')
    plt.title('Pathways Enriched in iPSC-derived vs. Primary Cardiomyocytes')
    plt.grid(axis='x')
    plt.tight_layout()
    

    # Save the plot to the desktop
    desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    plot_path = os.path.join(desktop_path, 'top_20_pathways.png')
    top_20_pathways_plot_path = plot_path
    plt.savefig(plot_path)

    # Extract the 'lead_genes' column from the first five rows of out_df
    lead_genes_first_five_rows = out_df.loc[:4, 'lead_genes']

    # Split the strings in 'lead_genes' column into lists of genes
    split_lead_genes = lead_genes_first_five_rows.str.split(';')

    # Create an empty set to store unique lead genes
    unique_lead_genes = set()

    # Iterate over each list of genes and add them to the set
    for genes_list in split_lead_genes:
        unique_lead_genes.update(genes_list)

    # Convert the set back to a sorted list
    all_lead_genes = sorted(list(unique_lead_genes))
    res['Gene'] = res.index
    volcano(res, symbol='Gene',log2fc_thresh = 2,colors = ['blue', 'dimgrey', 'red'])
    # Save the volcano plot to the desktop
    desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    volcano_plot_path = os.path.join(desktop_path, 'volcano_plot.png')
    plt.savefig(volcano_plot_path)

    # Define the layout for the PySimpleGUI window
    layout = [
        [sg.Text('Top 20 Pathways')],
        [sg.Image(key='-TOP_20_PATHWAYS-')],
        [sg.Text('Volcano Plot')],
        [sg.Image(key='-VOLCANO_PLOT-')],
        [sg.Button('OK')]
    ]

    # Create the PySimpleGUI window
    window = sg.Window('Top 20 Pathways and Volcano Plot', layout)

    window['-TOP_20_PATHWAYS-'].update(filename=top_20_pathways_plot_path)
    window['-VOLCANO_PLOT-'].update(filename=volcano_plot_path)
    
    # Event loop
    while True:
        event, values = window.read()
        if event == sg.WINDOW_CLOSED or event == 'OK':
            break

    # Close the window
    window.close()
    
    return log_output


if __name__ == '__main__':
    main()
