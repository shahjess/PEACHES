# PEACHES: Profiling Expression Analysis of Cellular Heterogeneity for Evaluating Stem Cell Differentiation (We can change this title, PEACHES is nice as an acronym for repo + environment title) 

This repo contains code for the PEACHES project (Team Members: Jessica Shah, Alana Mermin-Bunnell, Tess Fallon). 

# Data 
All data was sourced from GEO. **FILL IN ACCESSION NUMBERS**


# Folder structure 

```data```: 
/data 
contains all data directly used to build PEACHES platform
DE_ indicates output of a DESeq2 run used to build the platform 
log_output_ indicates output of a DESeq2 run used to test the platform 
essential_ and nonessential_ are files outputted by /code/get_essential_genes.py 

/data/raw
contains all raw data downloaded from GEO, categorized by sample type
this data was the data inputted into the DESeq2 for differential expression analysis 
 
```figures```: 
/Figures
contains .png files for all figures generated for the project

/figure_generation
contains .ipynb notebooks used to generate most figures in the project
we note that score_barGraph.png was generated from main() of diffScore.py

```code```: 
/code 
contains Python scripts used to develop and test the PEACHES scoring framework
contains script for running the GUI for the PEACHES framework 

# Installation and figure generation 

1. Clone the repository by running the following command in terminal:
```
git clone https://github.com/shahjess/PEACHES.git
```
2. Navigate current directory to the cloned repository directory in terminal:
```
cd PEACHES
```
3. To install the required Python modules, create and activate a [Conda](https://docs.anaconda.com/free/anaconda/install/) Environment 'PEACHES' by running the following commands in terminal: 
```
conda create -n PEACHES python=3.11  
conda activate PEACHES
conda install pip
pip install -r requirements.txt
```
4. 
