# SRUFinder

Detect and subtype Single Repeat Units, mini-arrays, and CRISPR arrays by repeat matching with a comprehensive repeat sequence database

# Table of contents
1. [Quick start](#quick)
2. [Installation](#install)
3. [SRUFinder](#sru)

## Quick start <a name="quick"></a>

```sh
conda create -n srufinder -c conda-forge -c bioconda -c russel88 srufinder
conda activate srufinder
srufinder my.fasta my_output
```

## Installation <a name="install"></a>
SRUFinder can be installed either through conda or pip.

It is advised to use conda, since this installs SRUFinder and all dependencies, and downloads with database in one go.

### Conda
Use [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/) to install.

Create the environment with SRUFinder and all dependencies and database
```sh
conda create -n srufinder -c conda-forge -c bioconda -c russel88 srufinder
```

### pip
If you have the dependencies (Python >= 3.8, Prodigal >= 2.6) in your PATH you can install with pip

Install srufinder python module
```sh
python -m pip install srufinder
```

Upgrade srufinder python module to the latest version
```sh
python -m pip install srufinder --upgrade
```


#### When installing with pip, you need to download the database manually: 
```sh
# Download
svn checkout https://github.com/Russel88/SRUFinder/trunk/data

# Tell SRUFinder where the data is:
# either by setting an environment variable (has to done for each terminal session, or added to .bashrc):
export SRUFinder_DB="/path/to/data/"
# or by using the --db argument each time you run SRUFinder:
srufinder input.fa output --db /path/to/data/
```

## SRUFinder - How to <a name="sru"></a>
SRUFinder takes as input a nucleotide fasta, and produces outputs with SRU predictions

#### Activate environment
```sh
conda activate srufinder
```

#### Run with a nucleotide fasta as input
```sh
srufinder genome.fa my_output
```

#### For metagenome assemblies and short contigs/plasmids/phages, change the prodigal mode
The default prodigal mode expects the input to be a single draft or complete genome
```sh
srufinder assembly.fa my_output --prodigal meta
```

#### Check the different options
```sh
srufinder -h
```

