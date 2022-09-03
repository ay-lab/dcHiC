# Re-selecting the PC

If a user doesn't like the selected chromosome pc for a sample he or she can manually force dchic to use their chosen pc using this code. 
This code has two options to reselect the pc. The first option 'ref' is if they have a separate signal file (mapped on same bedgraph as that of the hic data),
then reselectpc.r will correlate the bedGraph signal with the PC values and use it to select the PCs (e.g. re-run). The second option is manual re-selection 
(inputting the chromosome #, sample, PC #, etc.). 

Try `Rscript reselectpc.r --help` for all options. Below are some examples: 

- Automatically re-selecting: `Rscript reselectpc.r --reselect ref --rfile signal.txt`
- Manually re-selecting: `Rscript reselectpc.r --reselect man --sample <folder>_pca --chr <chr name> --pc <pc number> --pctype <intra or inter>`

# Using Existing PC values with dcHiC

Using `getcHiCinputfromExistingPCs.r`, anyone can reformat their existing pc values such that they can use dcHiC. The input must be a five column text file with the following fields: 
```
<pc bedGraph path>	<pc type>	<replicate name>	<sample name>
```

The first column is the path to the bedGraph file with PC information (4th column of the bedGraph file), the second column can have either 'intra' or 'inter' 
(whether component calculations were done with intra-chromosomal or inter-chromosomal data). The third column is the replicate name, and the fourth column
is the sample name. Replicate names and sample names should be different. 

# Time-series cluster the differential compartments
Using `pcacluster.r`, anyone can perform a time-series clustering of the normalized pc values from each sample. Check the help description for more details:
```
Options:
        --pcafile=PCAFILE
                dcHiC output file for example the differential.<intra/inter>_sample_group.Filtered.pcQnm.bedGraph file


        --samplefile=SAMPLEFILE
                A single column file that defines which columns to select from the bedGraph file


        --kcenter=KCENTER
                The number of centers to find from bedGraph file [Default 6]


        --genebed=GENEBED
                A sorted bed file with gene names in its fourth column


        --minprob=MINPROB
                Minimum probability for cluster membership [Default 0]


        --output=OUTPUT
                The prefix for all the output files [Default pc_clustered]


        ##### Note #####
        Please cite Wu M, Gu L (2021). TCseq: Time course sequencing data analysis. R package version 1.18.0 if you're using this code along with dcHiC


        -h, --help
                Show this help message and exit
```                

# Merge continuous differential compartments 

The mergeDEcompartments.sh script will help the user to merge continuous differential compartments. Use the following command to merge the differential compartments (FDR < 0.1 by default, change the value inside the script)

``` 
./mergeDEcompartments.sh -f <path/to/differential.intra_sample_group.pcOri.bedGraph> -d <distance in base pairs>
```
