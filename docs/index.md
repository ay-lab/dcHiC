# dcHiC HTML Demo Page

The following examples were created using dcHiC’s standalone visualization utility, which allows for single-command creation of HTML files for interactive compartment/genomic visualization. Examples from the paper are included below. Each file has a set of controls on the top, which allow users to search by region or by gene. These pages, available for download in the GitHub and created using IGV, will take a few seconds to load.

## ESC vs NPC (100kb)

See [this page](https://ay-lab.github.io/dcHiC/ESC_NPC_100Kb.RefineY_Rconf90_FDR10.pcOri.html) for a comparison of ES and NPC compartments at 100kb. It contains the following tracks:  
```markdown
ES/NPC compartments
Differential Compartment Scores
Differential Loops
LaminB1 Data
RNA-Seq Data (TPM)
HOMER Differential Scores (for comparison)
```

## ESC vs NPC (10kb)

See [this page](https://ay-lab.github.io/dcHiC/ESC_NPC_10Kb.RefineY_Rconf90_FDR10.pcOri.html) for a comparison of ES and NPC compartments at 10kb. It contains compartments and differential compartment scores. 

## ESC vs NPC vs CN (100kb)

See [this page](https://ay-lab.github.io/dcHiC/ESC_NPC_CN_100Kb.RefineY_Rconf90_FDR10.pcOri.html) for a complete mice neural differentiation overview. It contains the following tracks: 
```
ES/NPC compartments
Differential Compartment Scores
Differential Loops
LaminB1 Data
RNA-Seq Data (TPM)
```

### Notable Regions
In all three ESC vs NPC examples, notable **mESC**-specific genes include: *Dppa2/4* and *Zpf42*. Notable **NPC/CN**-specific genes include: *Ptn* and *Ephb1*. 

## Hematopoiesis 

See [this page](https://ay-lab.github.io/dcHiC/hematopoiesis_100Kb.RefineY_Rconf90_FDR10.pcOri.html) for compartment dynamics during mice hematopoetic differentiation. This page contains compartments and differential compartment scores. 

Notable changes include **Meis1**, **Runx2**, **Sox6**, and **Abca13**.

## Lympoblastoid Data

See [this page](https://ay-lab.github.io/dcHiC/gorking_hg19_40Kb.RefineY_Rconf90_FDR10.pcOri.html) for compartment scores across 20 lymphoblastoid cell lines, using data from [*Gorkin et al, 2019*](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1855-4). This page contains compartments and differential compartment scores. 

Notable differential regions include **NR2F2** and **PTPRK**.
