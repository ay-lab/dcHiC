# dcHiC: Differential Compartment Analysis of Hi-C Datasets [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

dcHiC is a tool for differential compartment analysis of Hi-C datasets. It employs Hierarchal Multiple Factor Analysis to normalize technical biases in two or more groups of Hi-C datasets within any hierarchal structure, before then using learned parameters from replicate data to call significant differential interactions in pairwise and group settings. Beyond this, dcHiC also has options to output HTML files for visualization (using IGV.js) and Gene Set Enrichment Analysis on significant compartment changes. 

For more information, see our presentation (Regulatory & Systems Genomics COSI) & <a href = "https://iscb-ismb20.myconferencenow.com/poster/dchic-differential-compartment-analysis-of-hi-c-datasets/"> poster</a> at ISMB 2020. 

## Installation

### Option 1: Conda (Highly Recommended)

We recommend using Conda to install all dependencies in a virtual environment. 

The suggested path is using the appropriate <a href="https://docs.conda.io/en/latest/miniconda.html/">Miniconda</a> distribution, which can be installed at the link. 

Make sure your "conda" command specifically calls the executable under the miniconda distribution (e.g., ~/miniconda3/condabin/conda). If "conda activate" command gives an error when you run it the first time then you will have to run "conda init bash" once.

To install, go to your directory of choice and run:

```bash
git clone https://github.com/ay-lab/dcHiC
conda env create -f ./dchic/environment.yml
conda activate dchic
```

Note: Conda has experienced some trouble as of late with Bioconductor packages. If Bioconductor-IHW raises an error, you can install all other packages and run the following command afterward within the activated virtual environment: 

```bash
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IHW")
```

### Option 2: Manual Installation

Clone the GitHub repository and install the following dependencies: 

Python Dependencies:
- python >= 3.6
- scipy
- numpy
- matplotlib
- scikit-learn

Other Dependencies: 
- R >= 3.4
- FactoMineR
- R Hashmap
- bedtools
- igv-reports 
- cooler (only if pre-processing _.cool_ files)

## Input

The input to dcHiC are tab-delimited O/E Hi-C correlation matrices. Pre-processing steps to convert files to this format (_.hic_, _.cool_, Hi-C Pro) are described later; however, whatever the format, you will have a directory structure as follows before entering to the application: 

```bash
./HiCDataset1
    chr1.matrix
    chr2.matrix
    ...
./HiCDataset2
    chr1.matrix
    chr2.matrix
    ...
...
```

dcHiC uses an input file with the format below. The replicate, name, and directory columns are required. The name column descibes a particular Hi-C profile, and each Hi-C profile must have two or more replicates (more replicates increases power of comparisons). Given these three columns, a two-tiered hierarchy will be run with all replicates under each "name" composing a group. Differential calling will be run between the average PC's of the replicates under each "name." 

```bash
replicate   name    (grouping)   directory
```

The "grouping" column is optional...


## TODO: Program Arguments, Examples, etc. 

