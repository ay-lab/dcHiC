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

Note: Conda has experienced some trouble as of late with Bioconductor packages. If Bioconductor-IHW raises an error during installation/while running, you may try one of several things. Linux is the recommended platform to perform analysis, but using a fresh miniconda installation can resolve library issues on Mac. Removing IHW from the yml file and install IHW directly via R may also work:

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
- igv-reports 

Other Dependencies: 
- R >= 3.4
- FactoMineR
- R Hashmap
- IHW (Bioconductor)
- bedtools
- Java 11+ (if you wish to perform gene set enrichment analysis)
- cooler (only if pre-processing _.cool_ files)
- R R.utils (if visualizing)

## About dcHiC

#### For detailed information on all of these methods, please see <a href = "https://www.dropbox.com/s/dpw2fcyx88un7y4/dcHiC%20Poster%20ISMB%20PPT%20FINAL.pdf?dl=0"> this poster </a> for specifics.

### Hierarchal Multiple Factor Analysis 

Multiple Factor Analysis is a variant of PCA (the traditional way to identify compartments) that normalizes biases between a cohort of datasets. dcHiC uses a hierarchal form of MFA that allows for comparison of any number of groups (in a two-tiered hierarchy) and datasets. 

### Differential Calling

Based on the groupings specified by the user in the input file (see below), dcHiC takes comparisons between Hi-C datasets in pairwise and group settings. By default, dcHiC takes every pairwise comparison and one "multi-comparison" across all groups. It outputs a "differential file" with all differential regions and a "full" file with all regions/values. 

### Visualization

Visualization is performed through IGV.js (package igv-reports), which creates standalone HTML files with built-in genome browsers. <a href = "https://dchic-viz.imfast.io/">See examples at this link</a>. For more information on how to run, see the input section below. 

### Compartmental Gene Set Analysis (cGSEA) 

To analyze the biological significance of differential compartments, dcHiC implements a pre-ranked GSEA (ranked by the -log10 of the p-values, m-distance, or dZsc). The gene set file should be specified by the user. This is usable but in development form. 

## Input File(s) Specifications

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

#### Important Note: We are aware it isn't ideal, but be sure names do not include underscores or hyphens (due to naming conventions in-program). Working on an update to resolve this issue. 

The optional "grouping" column can be thought of as an extra layer of organization. If you choose to include it, the same HMFA calculation will be run as before although dcHiC will take the average of all replicate PC values under each "grouping" rather than each "name" (which then averages different Hi-C profiles). <a href = "https://www.dropbox.com/sh/2lnsu3wz8j0gfz3/AAAG29_olvkRXuBcU4eFjJiTa?dl=0"> See sample input files here </a>. 

### Compartmental Gene Set Enrichment Analysis (cGSEA) Input File

cGSEA is an optional part of the pipeline that performs GSEA on a pre-ranked list of genes in significant differential compartments. To perform cGSEA, first ensure you have Java 11+. The input for cGSEA is a file
- The first line specifies how the data should be ranked: 1 is the -log(10) of the p-adjusted value, 2 is dZsc, and 3 is mahalanobis distance. In almost all cases, it should be 1 but <a href = "https://www.dropbox.com/s/dpw2fcyx88un7y4/dcHiC%20Poster%20ISMB%20PPT%20FINAL.pdf?dl=0"> this poster </a> provides specifics. 
- The second line is the path to a gene set GMT file. See the tutorial for an example. 
- The third line is the path to a gene position bed for the genome used. See <a href = "https://www.dropbox.com/sh/b8arrl7tzl1nzc3/AABV-1kSy93dB32peZ0ocfd4a?dl=0"> here </a> for examples.

## Program Arguments

To run dcHiC from top to bottom, use these arguments in dchic.py:

| Argument              | Meaning                 
| --------------------- | ----------------------- |
| **-res**                | Resolution of processing (i,e. 100000)
| **-inputFile**                | Path to input.txt file, as described above
| **-chrFile**                | File for chromosomes to be processed (one chromosome label per line) 
| **-input**                | Assign 1 (if using HOMER input) and 2 (for all else): Used to scan file names in input directories. 
| **-parallel**               | Optional: If you wish to use parallel processing for chromosomes, specify this option with the # of threads to be used. Otherwise, processing will be sequential by chromosome. 
| **-genome**         | Genome desired (hg38, hg19, mm10, mm9)
| **-signAnalysis**             | Specify which biological data will be used to determine eigenvector sign with ("gc" or "tss"). 
| **-alignData**             | Specify absolute path to UCSC goldenPath data to specify eigenvector sign. See <a href = "https://www.dropbox.com/sh/b9fh8mvkgbcugee/AABfzDQcF_Lt27TjfgrPswrta?dl=0">here</a> for examples. 
| **-cGSEA**             | Optional: If you wish to perform a ranked GSEA on signficant compartment changes, enter the input file as specified above.  
| **-keepIntermediates**             | Logical. Whether to keep certain intermediate files (such as R workspace data). Enter any argument (i,e. "1") to set true.

For instance, the following command would run a human non-HOMER input, 6 chromosomes at a time, with cGSEA, at a resolution of 100kb: 

```bash
python dchic.py -res 100000 -inputFile input.txt -chrFile chr.txt -input 2 -parallel 6 -genome hg38 -signAnalysis gc -alignData /path/to/hg38_goldenPathData -cGSEA cgsea.txt -keepIntermediates 1
```
This command would run a mice HOMER input, in sequence, without cGSEA or retaining intermediates, at a resolution of 500kb: 

```bash
python dchic.py -res 500000 -inputFile input.txt -chrFile chr.txt -input 1 -genome mm10 -signAnalysis gc -alignData /path/to/mm10_goldenPathData
```
dcHiC can be run from top to bottom or it can be run in a "modular" setting. See examples in the wiki (to be added).

## Visualization Input

Visualization can be run afterward using igvtrack.R, which takes two arguments. The first is a specified genome; the second is the input file. The first column should contain bedGraph files, the second column should define the name for each file, and the last column should contain labels for the group of each data. Visualization of compartment results _must_ use the "full_compartment_details" files in the DifferentialCompartment folder with the group name "compartment." 

```bash
file                     name          group
/path/to/full_bedGraph   expA_vs_expB  compartment
/path/to/full_bedGraph   expA_vs_expB  compartment
/path/to/other_bedGraph  laminB_expA   laminB
... 
```

With this, run:
```bash
Rscript /path/to/igvtrack.R [genome] [visualization file]
Rscript /path/to/igvtrack.R  mm10 viz.txt # an example
```

## Run dcHiC In Modular Fashion

Some run cases require a modular setup. For instance, if you want to run a large number (several dozens or more) of samples, it will be time-inefficient to run these in sequence with dchic.py and memory-intensive to run these in parallel. For these cases, dcHiC can be run in a modular fashion (by chromosome), a feature useful in environments like computing clusters where one simply needs to submit a job per chromosome. 

To run it for one particular chromosome XX, create a directory for that chromosome named "chr_XX" and use runhmfa.py (many of the same arguments as above). Rather than specifying the directory of input files in the last column, now specify the O/E correlation matrix file itself for that chromosome ("chrXX.matrix", for instance). 

After running several chromosomes separately, you can then go to the global directory that encompasses those directories and run differentialCalling.py (logical arguments), cgsea.py, and/or igvtrack.R. 

## Output

In the directory where you have run the data, you should have the following:

### A directory for each chromosome
Inside, the important files are: 
- O/E correlation matrices of common bins across input files (for experiment name XX):
```bash
BalancedChrMatrix_exp_XX.txt
```
- HMFA text and bedGraph results: "X" denotes experiment name, "XX" denotes experiment number in input, "XXX" denotes chromosome number
```bash
hmfa_X_exp_XX.txt
HMFA_chrXXX_exp_XX.bedGraph
```
- pcFiles directory: A directory of all raw PC files
- Other assorted files for program use

### A DifferentialCompartment directory
- All pairwise comparisons between Hi-C groups/groupings (specified in input) XXX and XXX:
```bash
    XXX_vs_XXX_full_compartment_details.bedGraph
    XXX_vs_XXX_differential_compartment_details.bedGraph
```
- Comparison file with significant differential interactions among all groups: 
```bash
    MultiComparison_full_compartment_details.bedGraph
    MultiComparison_differential_compartments_details.bedGraph
```
    
#### If cGSEA is to be performed
- Compartment detail files with genes outlined:
```
    Genes.XXX_vs_XXX_full_compartment_details.bedGraph
    ...
    Genes.MultiComparison_full_compartment_details.bedGraph
```
- Ranked gene files for input to GSEA:
```
    XXXvsXXX_genes_ranked.rnk
    MultiComparison_genes_ranked.rnk
```
- GSEA directories: Within each GSEA directory, there will be many results. View the index.html to explore results. See a guide to interpret results <a href = "https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm"> here </a>. 
```
    GSEA_XXX_vs_XXX
    ...
    GSEA_MultiComparison
```
  
### Other Information
- **Coordinate PNG's**: PC1 vs PC2 plots for experiment groups/groupings. 
- **PC Selection Info**: See which PC (1 or 2) was used for compartment analysis for each chromosome ("chr_info.txt")
            
## Tutorial: Mice Neural Differentiation Data

See <a href = "https://github.com/ay-lab/dcHiC/wiki/Mice-Neural-Differentiation-Tutorial">the tutorial page here</a>. 

## Contact

For help with installation or technical issues, contact:

Jeffrey Wang (jeffreywang@lji.org)

For general questions about the tool, interpretation, or technical details, contact:

Abhijit Chakraborty (abhijit@lji.org), Jeffrey Wang (jeffreywang@lji.org), Ferhat Ay (ferhatay@lji.org)
