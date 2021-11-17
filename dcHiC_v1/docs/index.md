# dcHiC HTML Demo Page

The following examples were created using dcHiC’s standalone visualization utility, which allows for single-command creation of beautiful, standalone HTML files that allow for interactive compartment/genomic visualization. Each file has a set of controls on the top, which allow users to search by region or by gene. These files, available for download in the GitHub and created using IGV Javascript, will take a few seconds to load.

See examples from the paper below, or build your own with this [tutorial](https://github.com/ay-lab/dcHiC/wiki/Mice-Neural-Differentiation-Tutorial).

## HOMER vs dcHiC

See [this page](https://ay-lab.github.io/dcHiC/html/dchic_homer.html) for a comparison of HOMER and dcHiC compartments. It contains the following tracks:  
```markdown
dcHiC Compartments
HOMER Compartments
mESC/NPC LaminB1 Profile
```

## Complete Mice Neural Differentiation (with dcHiC)

See [this page](https://ay-lab.github.io/dcHiC/html/multiWayMiceComplete.html) for complete mice neural differentiation overview. It contains the following tracks: 
```markdown
dcHiC Compartments
mESC/NPC LaminB1 Profile
RNA-Seq Data
ChIP Signals - H3K9me3, H3K4me3, H3K4me1, H3K27ac
```

Notable **mESC**-specific genes include: Dppa2/4, Zpf42   
Notable **NPC/CN**-specific genes include: Ptn, Ephb1, Dcx, Rmst   

## Hematopoiesis 

See [this page](https://ay-lab.github.io/dcHiC/html/hematopoiesis.html) for compartment dynamics during mice hematopoetic differentiation. This file also contains TPM values for 8 of the cell lines. 

Two notable genes with significant decreasing A-compartment changes include **Myc/Pvt1** and **Meis1**. Others include **Runx2, Sox6, and Abca13**.

## Human Tissue 

See [this page](https://ay-lab.github.io/dcHiC/html/NormalCellLines.html) for compartments across 11 normal human cell lines/tissue Hi-C samples. 

Notable **immune** (Naive Germinal B / GM12878) specific gene: CD84 
Notable **hippocampus**-specific gene: KLK6
