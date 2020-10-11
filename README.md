# dcHiC: Differential Compartment Analysis of Hi-C Datasets [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

dcHiC is a tool for differential compartment analysis of Hi-C datasets. It employs Hierarchal Multiple Factor Analysis to normalize technical biases in two or more groups of Hi-C datasets within any hierarchal structure, before then using learned parameters from replicate data to call significant differential interactions in pairwise and group settings. Beyond this, dcHiC also has options to output beautiful, standalone HTML files for visualization (using IGV.js) and several other useful analysis options. It is one of few tools that normalizes technical biases in Hi-C datasets and, to our knowledge, the only that performs Hi-C comparisons in group settings. 

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

### Filtering Structural Variations 

Large scale structural variations, like translocations and chromothripsis, create nonsense compartment calls. We introduce a measure "SVscore" to quantify this, and employ a filter as an option in dcHiC to address it. See the <a href = "https://github.com/ay-lab/dcHiC/wiki/Filtering-Chromosomes-With-Large-Structural-Variations">wiki page</a> for more.

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

#### Important Note: Be sure names do not include underscores or hyphens. 

The optional "grouping" column can be thought of as an extra layer of organization. If you choose to include it, the same HMFA calculation will be run as before although dcHiC will take the average of all replicate PC values under each "grouping" rather than each "name" (which then averages different Hi-C profiles). <a href = "https://www.dropbox.com/sh/2lnsu3wz8j0gfz3/AAAG29_olvkRXuBcU4eFjJiTa?dl=0"> See sample input files here </a>. 

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
| **-alignData**             | Specify absolute path to UCSC goldenPath data to specify eigenvector sign. See <a href = "https://www.dropbox.com/sh/b9fh8mvkgbcugee/AABfzDQcF_Lt27TjfgrPswrta?dl=0">here</a> for examples. If not included, dcHiC automatically downloads the necessary files. 
| **-keepIntermediates**          | Logical. Whether to keep certain intermediate files (such as R workspace data). Enter any argument (i,e. "1") to set true.
| **-blacklist**     |  Optional but HIGHLY recommended. Removes >1mb regions from the ENCODE blacklist before performing calculations. See "files" for hg19/hg38/mm9/mm10 blacklists. 
| **-ncp**   | The number of PC's to calculate & choose the final result from. Default is 2. Specify if more wanted.
| **-SVfilter** | Optional: If you wish to filter for structural variations, use the <a href = "https://github.com/ay-lab/dcHiC/wiki/Filtering-Chromosomes-With-Large-Structural-Variations">SVscore output</a> here. 
| **-removeFile** | Optional: If you do not wish to use the SVfilter option but simply want to remove a few chromosomes from a few samples, use this. Enter a tab-delimited file with experiment names on the left column (matching -inputFile) and chromosome numbers ('X', '2') on the right. 
| **-repParams** | Optional: If using data with no (or few) replicates, use a different set of replicate parameters instead. A set of high-quality parameters is available in the files for mice/human datasets. 

For instance, the following command would run a human non-HOMER input, with a blacklist, 6 chromosomes at a time, at a resolution of 100kb: 

```bash
python dchic.py -res 100000 -inputFile input.txt -chrFile chr.txt -input 2 -parallel 6 -genome hg38 -alignData /path/to/hg38_goldenPathData -keepIntermediates 1 -blacklist hg38blacklist_sorted.bed 
```
This command would run a mice HOMER input, in sequence, with a special replicate paramter file, without retaining intermediates, with a blacklist and SVfiltering, at a resolution of 500kb: 

```bash
python dchic.py -res 500000 -inputFile input.txt -chrFile chr.txt -input 1 -genome mm10 -alignData /path/to/mm10_goldenPathData -repParams mice_params.txt -SVfilter mice.svscore.txt -blacklist mm10blacklist_sorted.bed 
```

## Special Specifications

The blacklists are taken from a comprehensive study of problematic regions in the genome, dubbed the ENCODE blacklists. See the study <a href = "https://www.nature.com/articles/s41598-019-45839-z">here</a> and the full blacklists from the Boyle Lab <a href= "https://github.com/Boyle-Lab/Blacklist/tree/master/lists">here</a>.

Differential calling uses a large amalgmation of p-values across chromosomes to increase power. If more than 1/5 of chromosomes (specified in -chrFile) have some type of removal, either from SV filtering or manual removal, differential calling will instead be done on a chromosome-by-chromosome level. If there are only _some_ removals (in up to 1/5 of chromosomes), differential calling will be done chromosome-wise for those affected and together for the rest. 

The sample replicate parameter files for human datasets were created using Tier 1 ENCODE GM12878 and HMEC datasets for human samples (with no SV's detected via SVscore). The mice replicate parameter files were created using high-quality neural differentiation data (the same as that in our <a href = "https://github.com/ay-lab/dcHiC/wiki/Mice-Neural-Differentiation-Tutorial">tutorial</a>) with no SV's.  

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
  
### Other Information
- **Coordinate PNG's**: PC1 vs PC2 plots for experiment groups/groupings. 
- **PC Selection Info**: See which PC (1 or 2) was used for compartment analysis for each chromosome ("chr_info.txt" and "PCselection.txt")
            
## Tutorial: Mice Neural Differentiation Data

See <a href = "https://github.com/ay-lab/dcHiC/wiki/Mice-Neural-Differentiation-Tutorial">the tutorial page here</a>. 

## Updates

10/10: Substantial revisions made. Improved differential calling/PC selection, updated routines, SV-filtering, manual removal of chromosomes 

## Contact

For help with installation or technical issues, contact:

Jeffrey Wang (jeffreywang@lji.org)

For general questions about the tool, interpretation, or technical details, contact:

Abhijit Chakraborty (abhijit@lji.org), Jeffrey Wang (jeffreywang@lji.org), Ferhat Ay (ferhatay@lji.org)
