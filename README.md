# dcHiC: Differential Compartment Analysis of Hi-C Datasets [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

dcHiC is a tool for differential compartment analysis of Hi-C datasets. It employs Hierarchal Multiple Factor Analysis to normalize technical biases in two or more groups of Hi-C datasets within any hierarchal structure, before then using learned parameters from replicate data to call significant differential interactions in pairwise and group settings. Beyond this, dcHiC also has options to output beautiful, standalone HTML files for visualization (using IGV.js) and Gene Set Enrichment Analysis on significant compartment changes. It is one of few tools that normalizes technical biases in Hi-C datasets and, to our knowledge, the only that performs Hi-C comparisons in group settings. 

For more information, see our presentation (Regulatory & Systems Genomics COSI) & poster at ISMB 2020. This tool was produced by Jeffrey Wang, Abhijit Chakraborty, and Ferhat Ay at the La Jolla Institute for Immunology. 

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

## Input Specifications

The input to dcHiC are tab-delimited O/E Hi-C correlation matrices. Learn how to pre-process data (_.hic_, _.cool_, Hi-C Pro) on <a href = "https://github.com/ay-lab/dcHiC/wiki/Pre-Processing-Correlation-Matrices">the wiki page here</a>. Whatever option you choose, you should have directory structure like this before entering compartment analysis: 

```bash
./HiCDataset1.Rep1
    chr1.matrix
    chr2.matrix
    ...
./HiCDataset1.Rep2
    chr1.matrix
    chr2.matrix
    ...
./HiCDataset2.Rep1
    chr1.matrix
    chr2.matrix
    ... 
...
```

Create a file called input.txt for dcHiC with the format below. The replicate, name, and directory columns are required. The name column descibes a particular Hi-C profile, and each Hi-C profile must have two or more replicates (more replicates increases power of comparisons). Differential calling will be run between the average of the PC's of the replicates under each "name." 

```bash
replicate   name    (grouping)   directory
```

#### Important Note: Be sure names do not include underscores (due to naming conventions in-program)

The optional "grouping" column can be thought of as an extra layer of organization. If you choose to include it, the same HMFA calculation will be run as before. However, one important thing does change. Before, dcHiC would take the average of all replicate PC's under each "name" and make comparisons between the unique names in the second column. With a grouping (that encompasses multiple names), the average of all PC's under each grouping will be taken and comparisons will be made between those. See sample input files <a href = "https://www.dropbox.com/sh/2lnsu3wz8j0gfz3/AAAG29_olvkRXuBcU4eFjJiTa?dl=0> here </a>. 
    
## About dcHiC

### Hierarchal Multiple Factor Analysis 

Multiple Factor Analysis is a variant of PCA (the traditional way to identify compartments) that normalizes biases between a cohort of datasets. In short, it does so by performing a PCA on each individual dataset, dividing it by its first eigenvalue, then performing a global PCA on the cohort. HMFA is a hierarchal extension that performs this process over groupings of datasets in a hierachy. dcHiC implements HMFA for Hi-C datasets to allow for comparison of any number of groups (in a two-tiered hierarchy) and datasets. 

### Differential Calling

Based on the groupings specified by the user in the input file (see below), dcHiC takes comparisons between the average PC values of replicates for a cell line, or groups of cell lines. Differential analysis can be performed in pairwise and group settings. 

In short, our method first uses Mahalanobis Distance, a multivariate extension of z-score, to assign p-values to each bin based on its relative distance to the center of its distribution. After finding the outliers, it then takes pairwise comparisons of replicate PC values and learns data variation of replicates and uses those parameters to create a secondary measure "dZsc" that defines the degree to which bins across pairs of ._different cell lines._ differ differ beyond technical variation. For group analysis, comparisons are taken between all pairs and the dZsc values are averaged for each bin. 

These values, p-values and dZsc, are used together in an Independent Hypothesis Weighting (dZsc as the covariate) to create the output. By default, dcHiC takes every pairwise comparison and one "multi-comparison" across all groups (or groupings, if specified)—it outputs a "differential file" with all differential regions and a "full" file with all regions/values. 

#### Visualization

Visualization is performed through IGV.js (package igv-reports), which creates standalone HTML files with built-in genome browsers. To see examples, see this link: https://dchic-viz.imfast.io/. For more information on how to run, see the input section below. 

### Compartmental Gene Set Analysis

To analyze the biological significance of differential compartments, dcHiC implements a pre-ranked GSEA (ranked by the -log10 of the p-values, m-distance, or dZsc). The gene set file should be specified by the user. 

## Program Arguments

To run dcHiC from top to bottom, use these arguments in dchic.py:

| Argument              | Meaning                 
| --------------------- | ----------------------- |
| **-res**              | Resolution of processing (i,e. 100000)
| **-inputFile**              | Path to input.txt file, as described above
| **-chrFile**              | File for chromosomes to be processed (one chromosome label per line) 
| **-input**              | Assign 1 (if using HOMER input) and 2 (for all else): Used to scan file names in input directories. 
| **-parallel**             | Optional: If you wish to use parallel processing for chromosomes, specify this option with the # of threads to be used. Otherwise, processing will be sequential by chromosome. 
| **-genome**       | Genome desired (hg38, hg19, mm10, mm9)
| **-signAnalysis**           | Specify which biological data will be used to determine eigenvector sign with ("gc" or "tss"). 
| **-alignData**           | Specify absolute path to UCSC goldenPath data to specify eigenvector sign. See <a href = "https://www.dropbox.com/sh/odz8ietjutipbg9/AADt6y518gHo7ftPCxR-dZ0_a?dl=0">here</a> for examples. 
| **-cGSEA**           | Optional: If you wish to perform a ranked GSEA on signficant compartment changes, enter an input file where the first line is a number for GSEA ranking (1 = pAdj, 2 = dZsc, 3 = mahalanobis distance) and second line is the path to a gene set GMT file
| **-keepIntermediates**           | Logical. Whether to keep certain intermediate files (such as R workspace data). Enter any argument. 

For instance, the following command would run a human non-HOMER input, 6 chromosomes at a time, with cGSEA, at a resolution of 100kb: 

```bash
python dchic.py -res 100000 -inputFile input.txt -chrFile chr.txt -input 2 -parallel 6 -genome hg38 -signAnalysis gc -alignData /path/to/hg38_goldenPathData -cGSEA cgsea.txt -keepIntermediates 1
```
This command would run a mice HOMER input, in sequence, without cGSEA or retaining intermediates, at a resolution of 500kb: 

```bash
python dchic.py -res 500000 -inputFile input.txt -chrFile chr.txt -input 1 -genome mm10 -signAnalysis gc -alignData /path/to/mm10_goldenPathData
```
dcHiC can be run from top to bottom or it can be run in a "modular" setting. See examples in the wiki (to be added).

## Visualization Parameters

To run visualization, specify two parameters to igvtrack.R. The first is an input file. The first column should contain bedGraph files (full compartment details or other data), the second column should define the name for each file, and the last column should contain labels for the group of each data. Compartment details ._must._ be named "compartment." 

```bash
file                     name          group
/path/to/full_bedGraph   expA_vs_expB  compartment
/path/to/full_bedGraph   expA_vs_expB  compartment
/path/to/other_bedGraph  laminB_expA   laminB
... 
```

## Contact

For help with installation or technical issues, contact:

Jeffrey Wang (jeffreywang@lji.org)

For general questions about the tool, interpretation, or technical details, contact:

Abhijit Chakraborty (abhijit@lji.org), Jeffrey Wang (jeffreywang@lji.org), Ferhat Ay (ferhatay@lji.org)
