# dcHiC: Differential Compartment Analysis of Hi-C Datasets [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

dcHiC is a tool for differential compartment analysis of Hi-C datasets. It employs Hierarchal Multiple Factor Analysis to normalize technical biases in two or more groups of Hi-C datasets within any hierarchal structure, before then using learned parameters from replicate data to call significant differential interactions in pairwise and group settings. Beyond this, dcHiC also has options to output beautiful, standalone HTML files for visualization (using IGV.js) and Gene Set Enrichment Analysis on significant compartment changes. It is one of few tools that normalizes technical biases in Hi-C datasets and, to our knowledge, the only that performs Hi-C comparisons in groups settings. 

For more information, see our presentation (Regulatory & Systems Genomics COSI) & <a href = "https://iscb-ismb20.myconferencenow.com/poster/dchic-differential-compartment-analysis-of-hi-c-datasets/"> poster</a> at ISMB 2020. This tool was produced by Jeffrey Wang, Abhijit Chakraborty, and Ferhat Ay at the La Jolla Institute for Immunology. 

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
- Java 11+ (if you wish to perform gene set enrichment analysis)
- cooler (only if pre-processing _.cool_ files)

## About dcHiC

(To be finished)

## Input

The input to dcHiC are tab-delimited O/E Hi-C correlation matrices. Learn how to pre-process data (_.hic_, _.cool_, Hi-C Pro) on <a href = "https://github.com/ay-lab/dcHiC/wiki/Pre-Processing-Correlation-Matrices">the wiki page here</a>. Whatever option you choose, you should have directory structure like this before entering compartment analysis: 

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

Create a file called input.txt for dcHiC with the format below. The replicate, name, and directory columns are required. The name column descibes a particular Hi-C profile, and each Hi-C profile must have two or more replicates (more replicates increases power of comparisons). Differential calling will be run between the average of the PC's of the replicates under each "name." 

```bash
replicate   name    (grouping)   directory
```

The optional "grouping" column can be thought of as an extra layer of organization. If you choose to include it, the same HMFA calculation will be run as before. However, one important thing does change. Before, dcHiC would take the average of all replicate PC's under each "name" and make comparisons between the unique names in the second column. With a grouping (that encompasses multiple names), the average of all PC's under each grouping will be taken and comparisons will be made between those. Here are two such examples. 

## Program Arguments

To run dcHiC from top to bottom, use these arguments in dchic.py:

| Argument              | Meaning                 
| --------------------- | ----------------------- |
| **-res**              | Resolution of processing (i,e. 100000)
| **-inputFile**              | Path to input.txt file, as described above
| **-chrFile**              | File for chromosomes to be processed (one chromosome label per line) 
| **-input**              | Assign 1 (if using HOMER input) and 2 (for all else): Used to scan file names in input directories. 
| **-parallel**             | Optional: If you wish to use parallel processing for chromosomes, specify this option with the # of threads to be used
| **-genome**       | Genome desired (hg38, hg19, mm10, mm9)
| **-signAnalysis**           | Specify which biological data will be used to determine eigenvector sign with ("gc" or "tss"). 
| **-alignData**           | Specify absolute path to UCSC goldenPath data to specify eigenvector sign. See <a href = "https://www.dropbox.com/sh/odz8ietjutipbg9/AADt6y518gHo7ftPCxR-dZ0_a?dl=0">here</a> for examples. 
| **-cGSEA**           | Optional: dcHiC 
| **-keepIntermediates**           | Logical. Whether to keep certain intermediate files (such as R workspace data). Enter any argument. 

For instance: 

```bash
python dchic.py -res 100000 -inputFile input.txt -chrFile chr.txt -input 2 -parallel 6 -genome hg38 -signAnalysis gc -alignData /path/to/hg38_goldenPathData -cGSEA cgsea.txt -keepIntermediates 1
```

dcHiC can be run from top to bottom or it can be run in a "modular" setting. See examples in the wiki. (To be finished).

