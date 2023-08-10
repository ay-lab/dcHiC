# dcHiC: Differential Compartment Analysis of Hi-C Datasets 

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/274827002.svg)](https://zenodo.org/badge/latestdoi/274827002)

dcHiC is a tool for differential compartment analysis of Hi-C datasets. It features many capabilities, including:

- Optimized PCA calculations (faster + capable of analysis up to 5kb resolution)
- Comprehensive identification of significant compartment changes between any number of cell lines (with replicates), including with pseudo-bulk single cell data
- Beautiful standalone HTML files for visualization of results
- Identification of differential loops anchored in significant differential compartments (using [Fit-Hi-C](https://github.com/ay-lab/fithic))
- And much more! 

### Paper

If you want to see examples of dcHiC in action or cite our tool, please see our [paper in Nature Communications](https://www.nature.com/articles/s41467-022-34626-6)! Web-hosted visualization examples of case scenarios in the paper [here](https://ay-lab.github.io/dcHiC). 

To see how to run dcHiC, read our docs and try our demo (below)! Information about data pre-processing and running single-cell data is available in the [wiki](https://github.com/ay-lab/dcHiC/wiki). 

### Demo

This README contains the key information you will need to use this application. However, some users may find a demo helpful—ours includes a script to run package installation as well as detailed guides for different options of dcHiC. All of these resources are available in `packages/dchic_demo.zip`, with relevant instructions inside! 

## Installation

The latest version of dcHiC runs pre-dominantly from R (3+) and Python (3+). The necessary packages may be installed via conda or manually (those transitioning environments should have most, if not all, of the packages already installed). For the core application, the following packages are necessary:

### Option 1: Conda

We recommend using Conda to install all dependencies in a virtual environment. The suggested path is using the appropriate <a href="https://docs.conda.io/en/latest/miniconda.html/">Miniconda</a> distribution. 

If you face any issues, be sure your "conda" command specifically calls the executable under the miniconda distribution (e.g., ~/miniconda3/condabin/conda). If "conda activate" command gives an error when you run it the first time then you will have to run "conda init bash" once.

To install, go to the directory of your choice and run:

```bash
git clone https://github.com/ay-lab/dcHiC
conda env create -f ./packages/dchic.yml
conda activate dchic
```

Afterward, activate the environment and install some purpose-built processing functions with `R CMD INSTALL functionsdchic_1.0.tar.gz` (functions file under 'packages'). M1 Mac users may face some issues, as some bioconductor packages have not yet been updated for native ARM64 support; we recommend using an x86-64 based OS for the cleanest experience. 

### Option 2: Manual Installation

To install the dependencies manually, ensure that you have the following packages installed:

### Packages in R
- Rcpp
- optparse
- bench
- bigstatsr
- bigreadr
- robust
- data.table
- networkD3
- depmixS4
- rjson
- limma (bioconductor)
- IHW (bioconductor)
- R.utils
- hashmap (.tar.gz file under 'packages')

### Packages in Python
- igv-reports

### Bedtools 
- dcHiC requires [bedtools](https://bedtools.readthedocs.io/en/latest/). Please install the program as directed—it should be accessible via $PATH. 

Those who wish to perform differential loop analysis should also download the latest Python version of FitHiC, which requires a set of [Python libraries](https://github.com/ay-lab/fithic): numpy, scipy, sk-learn, sortedcontainers, and matplotlib. You may also need to install 'cooler' if you wish to use *.cool* files. See [documentation](https://github.com/ay-lab/dcHiC/wiki/Pre-Processing-Data) on how to do so. 

Afterward, activate the environment and install some purpose-built processing functions with `R CMD INSTALL functionsdchic_1.0.tar.gz` (functions file under 'packages'). 

### To check which R packages are already installed
```
Rscript -e 'plist <- c("functionsdchic","hashmap","R.utils","Rcpp","RcppEigen","BH","optparse","bench","bigstatsr","bigreadr","robust","data.table","networkD3","depmixS4","rjson","limma","IHW"); setdiff(plist,basename(find.package(plist)))'
````
If you get character(0) then you're all set, otherwise install the packages shown in the output.

## Input File

Create an input file for dcHiC with the format below. The matrix and bed columns are for input data (see next section), whereas the replicate_prefix and experiment_prefix columns describe the hierarchy of data. 

**Note: Do not use dashes ("-") or dots (".") in the replicate or experiment prefix names.**

```bash
<mat>         <bed>         <replicate_prefix>      <experiment_prefix>
```

For instance, consider this sample file which describes two replicates for two Hi-C profiles:

```
matr1_e1.txt  matr1_e1.bed   exp1_R1_100kb                  exp1
matr2_e1.txt  matr2_e2.bed   exp1_R2_100kb                  exp1
matr1_e2.txt  matr1_e2.bed   exp2_R1_100kb                  exp2
matr2_e2.txt  matr2_e2.bed   exp2_R2_100kb                  exp2
```

## Input Data

dcHiC accepts sparse matrices as its input (Hi-C Pro style). If you have *.cool* or *.hic* files, see how to convert their format [here](https://github.com/ay-lab/dcHiC/wiki/Pre-Processing-Data). 

To see the full list of options, run `Rscript dchicf.r --help` or view `dchicdoc.txt` [here](https://github.com/ay-lab/dcHiC/blob/master/dchicdoc.txt).

The matrix file should look like this:

```
<indexA> <indexB> <count>

1         1       300
1         2       30
1         3       10
2         2       200
2         3       20
3         3       200
 			....
```

... And the corresponding bed file like this:

```
<chr>	<start>	<end>	<index>

chr1	0	      40000	   1
chr1	40000	  80000	   2
chr1	80000	  120000	 3
 			....
```

### Blacklisted Regions

Many high-throughput genomics studies "blacklist" problematic mapping regions (see the study <a href = "https://www.nature.com/articles/s41598-019-45839-z">here</a>). If you wish to blacklist regions from your data, you may do so by adding a fifth column to your input file containing 1's in rows that should be blacklisted:

```
<chr>	<start>	<end>	<index>	<blacklisted>

chr1	0	      40000	 1	     0
chr1	40000	  80000	 2	     1
 			....
```

## Run Options

To see the full list of run options with examples of run code for each one, run `Rscript dchicf.r --help`. The most high-level option is `--pcatype`, which allows users to perform different types of step-wise analysis. Each of these run options will require other input information.

| --pcatype option              | Meaning                 
| --------------------- | ----------------------- |
| **cis**                | Find compartments on a cis interaction matrix
| **trans**                | Find compartments on a trans interaction matrix
| **select**                | Selection of best PC for downstream analysis [Must be after cis or trans step]
| **analyze**                | Perform differential analysis on selected PC's [Must be after select step]
| **subcomp**                | Optional: Assigning sub-compartments based on PC magnitude values using HMM segmentation 
| **fithic**         | Run [Fit-Hi-C](https://github.com/ay-lab/fithic) to identify loops before running dloop (Optional but recommended)
| **dloop**             | Find differential loops anchored in at least one of the differential compartments across the samples (Optional but recommended)
| **viz**  | Generate IGV vizualization HTML file. Must have performed other steps in order (optional ones not strictly necessary) before this one.
| **enrich**     |  Perform gene enrichment analysis (GSEA) of genes in differential compartments/loops


Here is a sample full run using the traditional cis matrix for compartment analysis: 

```
Rscript dchicf.r --file input.ES_NPC.txt --pcatype cis --dirovwt T --cthread 2 --pthread 4
Rscript dchicf.r --file input.ES_NPC.txt --pcatype select --dirovwt T --genome mm10
Rscript dchicf.r --file input.ES_NPC.txt --pcatype analyze --dirovwt T --diffdir ES_vs_NPC_100Kb
Rscript dchicf.r --file input.ES_NPC.txt --pcatype fithic --dirovwt T --diffdir ES_vs_NPC_100Kb --fithicpath "/path/to/fithic.py" --pythonpath "/path/to/python"
Rscript dchicf.r --file input.ES_NPC.txt --pcatype dloop --dirovwt T --diffdir ES_vs_NPC_100Kb
Rscript dchicf.r --file input.ES_NPC.txt --pcatype subcomp --dirovwt T --diffdir ES_vs_NPC_100Kb
Rscript dchicf.r --file input.ES_NPC.txt --pcatype viz --diffdir ES_vs_NPC_100Kb --genome mm10 
Rscript dchicf.r --file input.txt --pcatype enrich --genome mm10 --diffdir conditionA_vs_conditionB --exclA F --region both --pcgroup pcQnm --interaction intra --pcscore F --compare F
```

## Output

As output, dcHiC creates two types of directories. The first are raw PCA results, in directories named after the third column of the input file. One of these is created for each input Hi-C profile; inside, there will be directories "intra_pca" or "inter_pca" depending on whether the user specified compartment calculations based on intra- or inter-chromosomal interactions and raw PC values for each chromosome inside each one. 

The second overarching directory is called `DifferentialResult`, which contains directories for differential results (on any number of parameter settings). These directory names are specified under the `-analyze` pcatype option (which performs differential calling) dcHiC where users denote a `--diffdir` where they want the analysis to be done. Multiple directories, with different analysis parameters, can be stored under the global DifferentialResult directory. 

Inside each *diffdir*, there are raw compartment results ("expXX_data") and two PC output directories `PcOri` and `PcQnm` with combined and quantile-normalized compartment results. Finally, there will be a directory `fdr_result` containing differential compartment, loop, and subcompartment results. Inside `fdr_result`, the `sample_combined` files contain complete bedGraphs with average PC values across replicates for all *XX* cell lines, as well as a final adjusted p-value denoting the significance of changes between Hi-C experiments for that compartment bin. The `sample_combined.Filtered` files contain the same information, filtered by a p-value cutoff. 

Other `subcompartments` and `compartmentLoops` may be there depending on whether the user opted to run those options. The differential loop files list significant loop interactions and their associated differential compartment anchors, whereas the `subcompartment` files illustrate HMM-segmented subcompartments based on the magnitude of the PC values. 

Below is a diagram of the overarching results structure, containing two different runs (
```
dcHiC_dir
 exp1_rep1_100kb_pca
   intra_pca
      [files]
   inter_pca
      [files]
 exp1_rep2_100kb_pca
 exp2_rep1_100kb_pca
 exp2_rep2_100kb_pca
 DifferentialResult
   inter_100kb_diff
     [files]
   intra_100kb_diff
     exp1_data
     exp2_data
     fdr_result
     fithic_run
     geneEnrichment
     pcOri
     pcQnm
     viz
```

## Technical Specifications / Errors To Check 

There are a few technical implementation items to note:

**Chromosomes:** If you are running into issues during running dcHiC, removing chrM, chrY and other non-standard chromosomes will help. There have been many issues raised about this; we highly recommend you [search](https://github.com/ay-lab/dcHiC/issues?q=is%3Aissue+label%3A%22user+questions%22+is%3Aclosed) for the label "user questions" under Issues if you encounter an error related to this. Also make sure that the chromosome labels in the matrices match the goldenPath files; see [this issue](https://github.com/ay-lab/dcHiC/issues/47). 

**`fithic`/`dloop`**: If running `dloop`, dcHiC will first run Fit-Hi-C on the data. You will need to follow the Fit-Hi-C running procedure to do this, which will require generating a bias file. See "FitHiC2 bias file format" [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7451401/). 

**Support for other genomes:** While it has only been extensively tested for human and mouse genomes, dcHiC supports most other commonly-used genomes that are under the UCSC [genome page](https://hgdownload.soe.ucsc.edu/goldenPath/). To utilize this, create a folder `*{genome}_{resolution}_goldenpathData*` (e.g hg38_100000_goldenpathData).

Within that folder put three files: 
- `{genome}.fa` (e.g. [hg38.fa](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)) 
- `{genome}.tss.bed` (e.g. [hg38.tss.bed](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz), the TSS file. Please make sure the TSS position is selected based on the strad direction!) Note that this may be named `.refGene.gtf.gz`. 
- `{genome}.chrom.sizes` (e.g. [hg38.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)). 

These files can be found under the UCSC bigZips page for the specified genome. When running dcHiC use the `--gfolder` option in the `select` step to provide the folder path, and dcHiC will create the necessary files. 
 
**Compartment clustering:** Due to statistical noise, edge cases, and other factors, lone differential compartments occassionally crop up (ex: one bin is "significant" but all of its neighbors are not). These may be significant if analyzing at coarse resolution, but can also be misleading, especially if analyzing at very fine resolution. By default, dcHiC does not filter any of these lone compartments; however, there are two parameters to do so: `distclust` is the distance threshold for close differential regions to be a "cluster." If it's 0, only adjacent differential compartments form a cluster. If it's 1, differential compartments separated by up to 1 bin are a cluster. The other parameter is `numberclust`, which is a filter for the minimum number of significant bins within a cluster. 

**Chromosome-arm wise PCA calculation**: In order to perform p and q-arm wise PCA calculations, please check the __run_dcHiC_chrArms_pca_step1.pl__ and __run_dcHiC_chrArms_combine_step2.pl__ scripts provided under the ```utility/Chromosome_ArmWise_PCA/``` folder.  

## dcHiC v1

We previously released a different version of dcHiC (under the branch "dcHiC-v1") based on Python & R. While we hope that all users try the latest version of dcHiC, all code and documentation for the first version remains and we will continue offering support for it into the future. 

## Contact

For help with installation, technical issues, interpretation, or other details, feel free to raise an issue or contact us: 

Abhijit Chakraborty (abhijit@lji.org), Jeffrey Wang (jeffreywang@lji.org), Ferhat Ay (ferhatay@lji.org)

## Citation 

[Chakraborty A, Wang JG, Ay F. dcHiC detects differential compartments across multiple Hi-C datasets. Nat Commun. 2022;13(1):6827. Published 2022 Nov 11. doi:10.1038/s41467-022-34626-6](https://www.nature.com/articles/s41467-022-34626-6)
