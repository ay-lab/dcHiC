# Chromosome-arm wise PCA calculation

Presence of centromeres in the human chromosomes (and others) may hinder the calculation compartment scores and thereby subsequent differential analysis. In order to avoid such scenario, a user may choose to calculate the p and q-arm compartment scores separately. The two following scripts will help to calculate the chromosome-arm wise compartment scores. 

The __run_dcHiC_chrArms_pca_step1.pl__ script will extract the p and q-arms for each chromosome invidually and create the required input files for dcHiC. Once the first step of dcHiC is finished, the second script __run_dcHiC_chrArms_combine_step2.pl__ will combine the arm wise compartment scores for each chromosome for the subsequent differential analysis by dcHiC. The user should provide an input.txt file and a centromeres.txt file to the __run_dcHiC_chrArms_pca_step1.pl__ script. As an example we are providing the human (hg38) centromeres.txt file. 

```
perl run_dcHiC_chrArms_pca_step1.pl

HELP:
        --cpu=Number of CPUs to be used for the job
        --mem=Amount of memory to be used for the job (in GB)
        --hrs=Amount of wall time requested for the job (in hours)
        --pfx=Prefix of job name
        --cmd=input.txt file
        --pexcl=Exclude chromosomes from parm pca analysis (provide multiple chromosome names as e.g. chr21,chr22)
        --qexcl=Exclude chromosomes from qarm pca analysis (provide multiple chromosome names as e.g. chr21,chr22)
        --cen=centromere file (Downlod centromere file from https://www.ncbi.nlm.nih.gov/grc/human)
```
```
perl run_dcHiC_chrArms_combine_step2.pl

HELP:
        --cmd=input.txt file
```
