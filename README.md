# attack_amr
Bioinformatics pipeline for the antimicrobial resistance (AMR) analysis of hospital wastewater metagenomic samples

## 0. Introduction
Welcome to ATTACK-AMR - a bioinformatics pipeline for analysis of hospital wastewater antimicrobial resistance metagenomic data.

## 1. Installation
This pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
Additional dependencies not handled by Snakemake are described in Section 1.3.

### 1.1. Install Miniconda 
```
$ curl -sL \
  "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > \
  "Miniconda3.sh"
$ bash Miniconda3.sh
$ conda update conda
$ rm Miniconda3.sh
$ conda install wget
```

### 1.2. Install Mamba 
```
$ conda config --add channels conda-forge
$ conda update -n base --all
$ conda install -n base mamba
```

### 1.3. Install Snakemake
```
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This creates an isolated enviroment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test snakemake:
```
$ snakemake --help
```

### 1.4. Install Additional Dependencies
Install git and gawk. We require gawk to process the filtering stage of our databases.
```
$ mamba install git
$ mamba install gawk
```

### 1.5. Download the pipeline
Download ATTACK-AMR from the online [repository](https://github.com/bioinfodlsu/attack_amr), or using the command line:
```
$ git clone https://github.com/bioinfodlsu/attack_amr
```

## 2. Quickstart Usage Guide

### 2.1. Input
The pipeline requires, at the very least: (1) Metagenomic sequences (sample sequences can be downloaded at [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB47975), and (2) reference protein databases for ([Resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/), [MGE](https://github.com/KatariinaParnanen/MobileGeneticElementDatabase)).
These and other input parameters are specified via a YAML-format config file -- config.yaml is provided in the config folder. 

### 2.2. Running the pipeline
After constructing a config.yaml file and with the snakemake conda environment you created earlier activated, you can call the pipeline from the top-level directory of ATTACK-AMR:
```
$ cd attack_amr 
$ snakemake --use-conda --cores all
```

### 2.3. Ouput
Outputs are stored the top-level directory of ATTACK-AMR.
