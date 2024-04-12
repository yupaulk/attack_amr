# attack_amr
Source code for analysis of hospital wastewater metagenomic data

# 0. Introduction
Welcome to ATTACK-AMR - a bioinformatics pipeline for analysis of hospital wastewater antimicrobial resistance metagenomic data.

# 1. Installation
This pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
Additional dependencies not handled by Snakemake are described in Section 1.3.

## 1.1. Install Conda 
Download Miniconda3 installer for Linux from  [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
Installation instructions are [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).
Once installation is complete, you can test your Miniconda installation by running:
```
$ conda list
```

## 1.2. Install Snakemake
Snakemake recommends installation via Conda:
```
$ conda install -c conda-forge mamba
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

## 1.3. Install Additional Dependencies
We require gawk to process the filtering stage of our databases.
```
$ mamba install gawk
```
To test gawk:
```
$ gawk --version
```

## 1.4. Download the pipeline
Download ATTACK-AMR from the online [repository](https://github.com/bioinfodlsu/attack_amr), or using the command line:
```
git clone https://github.com/bioinfodlsu/attack_amr
```

# 2. Quickstart Usage Guide

# 2.1. Input
The pipeline requires, at the very least: (1) Metagenomic sequences (sample sequences can be downloaded at [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB47975), and (2) reference protein databases for ([Resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/), [MGE](https://github.com/KatariinaParnanen/MobileGeneticElementDatabase)).
These and other input parameters are specified via a YAML-format config file -- config.yaml is provided in the config folder. 

# 2.2. Running the pipeline
After constructing a config.yaml file and with the snakemake conda environment you created earlier activated, you can call the pipeline from the top-level directory of ATTACK-AMR:
```
cd attack_amr 
snakemake --use-conda --cores all
```

# 2.3. Ouput
Outputs are stored the top-level directory of ATTACK-AMR.
