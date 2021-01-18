# dcHiC: Differential Compartment Analysis of Hi-C Datasets [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

dcHiC is a tool for differential compartment analysis of Hi-C datasets. It employs Multiple Factor Analysis to normalize technical biases in two or more groups of Hi-C datasets within any hierarchal structure, before then using learned parameters from replicate data to call significant differential interactions in pairwise and group settings. Beyond this, dcHiC also has options to output beautiful, standalone HTML files for visualization and several other useful analysis options (web-hosted examples can be seen [here](https://ay-lab.github.io/dcHiC)). It is one of few tools that normalizes technical biases in Hi-C datasets and, to our knowledge, the only that performs Hi-C comparisons in group settings. 

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
- R R.utils
- cooler (only if pre-processing _.cool_ files)

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

The chromosome text file should be one chromosome label per line like follows:
```bash
1
2
3
...
X
```

Create a file called input.txt for dcHiC with the format below. The replicate, name, and directory columns are required. The name column descibes a particular Hi-C profile, and each Hi-C profile must have two or more replicates (more replicates increases power of comparisons). Differential calling will be run between the average of the PC's of the replicates under each "name." 

```bash
replicate   name    (grouping)   directory
```

For instance:
```bash
GM12878.0   GM12878 /path/to/GM_0
GM12878.1   GM12878 /path/to/GM_1
GM12878.2   GM12878 /path/to/GM_2
HMEC.0  HMEC    /path/to/HMEC_0
HMEC.1  HMEC    /path/to/HMEC_1
```

#### Important Note: Be sure names do not include underscores or hyphens. 

The optional "grouping" column can be thought of as an extra layer of organization. If you choose to include it, the same MFA calculation will be run as before although dcHiC will take the average of all replicate PC values under each "grouping" rather than each "name" (which then averages different Hi-C profiles). <a href = "https://www.dropbox.com/sh/2lnsu3wz8j0gfz3/AAAG29_olvkRXuBcU4eFjJiTa?dl=0"> See sample input files here </a>. 

## Program Arguments

To run dcHiC from top to bottom, use these arguments in dchic.py:

| Argument              | Meaning                 
| --------------------- | ----------------------- |
| **-res**                | Resolution of processing (i,e. 100000)
| **-inputFile**                | Path to input.txt file, as described above
| **-chrFile**                | File for chromosomes to be processed (one chromosome label per line)
| **-input**                | Assign 1 (if using HOMER input) and 2 (for all else)
| **-parallel**               | Optional: If you wish to use parallel processing for chromosomes, specify this option with the # of threads to be used. Otherwise, processing will be sequential by chromosome. 
| **-genome**         | Genome desired (hg38, hg19, mm10, mm9)
| **-alignData**             | Specify absolute path to UCSC goldenPath data to specify eigenvector sign. See <a href = "https://www.dropbox.com/sh/b9fh8mvkgbcugee/AABfzDQcF_Lt27TjfgrPswrta?dl=0">here</a> for examples. If not included, dcHiC automatically downloads the necessary files. 
| **-keepIntermediates**  | Logical. Whether to keep certain intermediate files (such as R workspace data). Enter any argument (i,e. "1") to set true.
| **-blacklist**     |  Optional but HIGHLY recommended. Removes >1mb regions from the ENCODE blacklist before performing calculations. See "files" for hg19/hg38/mm9/mm10 blacklists. 
| **-ncp**   | The number of PC's to calculate & choose the final result from. Default is 2. Specify if more wanted.
| **-repParams** | Optional: If using data with no (or few) replicates, use a different set of replicate parameters instead. A set of high-quality parameters is available in the files for mice/human datasets. 

For instance, the following command would run a human non-HOMER input, with a blacklist, 6 chromosomes at a time, at a resolution of 100kb: 

```bash
python dchic.py -res 100000 -inputFile input.txt -chrFile chr.txt -input 2 -parallel 6 -genome hg38 -alignData /path/to/hg38_goldenPathData -keepIntermediates 1 -blacklist hg38blacklist_sorted.bed 
```
This command would run a mice HOMER input, with chromosomes in sequence, without retaining intermediates, with a blacklist, at a resolution of 500kb: 

```bash
python dchic.py -res 500000 -inputFile input.txt -chrFile chr.txt -input 1 -genome mm10 -alignData /path/to/mm10_goldenPathData -blacklist mm10blacklist_sorted.bed 
```

## Special Specifications

### Running dcHiC Without Replicates

Differential calling with dcHiC "learns" the amount that PC (compartment) values vary between biological replicate datasets and uses those parameters for significance thresholds. However, it is also possible to run dcHiC without replicates. Using the run option "-repParams," one can use a pre-trained parameter file (in "files") like this: 

```bash
python dchic.py -res 100000 -inputFile input.txt -chrFile chr.txt -input 2 -parallel 4 -genome mm10 -alignData /path/to/mm10_goldenPathData -repParams miceparams.txt -blacklist mm10blacklist_sorted.bed 
```

The sample replicate parameter files for human datasets were created using Tier 1 ENCODE GM12878 and HMEC datasets for human samples/ The mice replicate parameter files were created using gold standard mice neural differentiation data (the same as that in our <a href = "https://github.com/ay-lab/dcHiC/wiki/Mice-Neural-Differentiation-Tutorial">tutorial</a>). 


### Standalone Differential Calling

Issues can arise post-processing (for instance: poor concordance in PC values between replicates). In these cases, it is possible to run the differential calling segment standalone as well. The arguments are straightforward and can be accessed with "-h".  Here is one example: 

```bash
python differentialCalling.py -inputFile input.txt -chrFile chr.txt -multiComp 1 -res 100000 -blacklist mm10blacklist_sorted.bed -genome mm10 -repParams miceparams.txt
```

### Other Notes 

Genome blacklisted regions are taken from a comprehensive study of problematic regions in high-throughput sequencing experiments, dubbed the ENCODE blacklists. These are available in the "files" directory. See the study <a href = "https://www.nature.com/articles/s41598-019-45839-z">here</a> and the full blacklists from the Boyle Lab <a href= "https://github.com/Boyle-Lab/Blacklist/tree/master/lists">here</a>. 

## Visualization Input

Visualization can be run afterward using the igvtrack R script, which takes two arguments. The first is a specified genome; the second is the input file. The first column should contain bedGraph files, the second column should define the name for each file, and the last column should contain labels for the group of each data. Visualization of compartment results _must_ include the "full_compartment_details" files in the DifferentialCompartment folder with the group name "compartment." 

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

See examples at [ay-lab.github.io/dcHiC](https://ay-lab.github.io/dcHiC). The dark green track (p-adjusted values) represent the FDR-adjusted p-values. 

## Output

In the directory where you have run the data, you should have the following:

### A directory for each chromosome
Inside, the important files are: 
- O/E correlation matrices of common bins across input files (for experiment name X):
```bash
BalancedChrMatrix_exp_X.txt
```
- Text and bedGraph results: "X" denotes experiment name, "XX" denotes experiment number, "XXX" denotes chromosome number
```bash
hmfa_X_exp_XX.txt
HMFA_chrXXX_exp_X.bedGraph
```
- A directory of all raw PC files
- Other assorted files for program use

### A Differential Compartment directory
- All pairwise comparisons between Hi-C groups/groupings (specified in input) XXX and XXX:
```bash
    XXX_vs_XXX_full_compartment_details.bedGraph
    XXX_vs_XXX_differential_compartment_details.bedGraph
```
- A multi-way comparison file with significant differential interactions among all groups: 
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

11/14: Bug fixes and slight routine tweaks. Updated blacklists. 

10/10: Substantial revisions made. Improved differential calling/PC selection, updated routines, SV-filtering, manual removal of chromosomes options. 

## Contact

For help with installation or technical issues, contact:

Jeffrey Wang (jeffreywang@lji.org)

For general questions about the tool, interpretation, or technical details, contact:

Jeffrey Wang (jeffreywang@lji.org), Abhijit Chakraborty (abhijit@lji.org), Ferhat Ay (ferhatay@lji.org)
