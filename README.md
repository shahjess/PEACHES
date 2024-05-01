# PEACHES: Profiling Expression Analysis of Cardiomyocyte Heterogeneity for Evaluating Stem Cell Differentiation 

This repo contains code for the PEACHES project (Team Members: Jessica Shah, Alana Mermin-Bunnell, Tess Fallon). 

# Data 
Publicly available human bulk RNA-seq datasets were curated from the NCBI GEO of iPSC-derived cardiomyocytes, ESC-derived cardiomyocytes, fibroblast-derived cardiomyocytes, undifferentiated iPSCs, and iPSC-derived lung epithelial cells (Table S1). Human heart left ventricle RNA-seq data was obtained from five samples (Sample ID: ENSG00000223972.5, ENSG00000227232.5, ENSG00000278267.1, ENSG00000243485.5, ENSG00000237613.2) in the Adult Genotype-Tissue Expression (GTEx) v8 database. Gene expression values from these samples were derived from non-diseased normal tissues of postmortem human donors (aged 21–70).

# Folder structure 

```/data```: 
contains all data directly used to build PEACHES platform
```DE_``` indicates output of a DESeq2 run used to build the platform 
```log_output_``` indicates output of a DESeq2 run used to test the platform 
```essential_``` and ```nonessential_``` are files outputted by ```/code/get_essential_genes.py``` 

```/data/raw```:
contains all raw data downloaded from GEO, categorized by sample type
this data was the data inputted into the DESeq2 for differential expression analysis 
 
```/figures```: 
contains .png files for all figures generated for the project

```/figure_generation```:
contains .ipynb notebooks used to generate most figures in the project.
Note that score_barGraph.png was generated from main() of diffScore.py

```/code```: 
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
