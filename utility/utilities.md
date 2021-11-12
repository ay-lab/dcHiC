# Re-selecting the PC

If a user doesn't like the selected chromosome pc for a sample he or she can manually force dchic to use their chosen pc using this code. 
This code has two options to reselect the pc. The first option 'ref' is if they have a separate signal file (mapped on same bedgraph as that of the hic data),
then reselectpc.r will correlate the bedGraph signal with the PC values and use it to select the PCs (e.g. re-run). The second option is manual re-selection 
(inputting the chromosome #, sample, PC #, etc.). 

Try `Rscript reselectpc.r --help` for all options. Below are some examples: 

- Automatically re-selecting: `Rscript reselectpc.r --reselect ref --rfile signal.txt`
- Manually re-selecting: `Rscript reselectpc.r --reselect man --sample <folder>_pca --chr <chr name> --pc <pc number> --pctype <intra or inter>`

# Using Existing PC values with dcHiC

Using this code, anyone can reformat their existing pc values such that they can use dcHiC. The input must be a five column text file with the following fields: 
```
<pc bedGraph path>	<pc type>	<replicate name>	<sample name>
```

The first column is the path to the bedGraph file with PC information (4th column of the bedGraph file), the second column can have either 'intra' or 'inter' 
(whether component calculations were done with intra-chromosomal or inter-chromosomal data). The third column is the replicate name, and the fourth column
is the sample name. Replicate names and sample names should be different. 
