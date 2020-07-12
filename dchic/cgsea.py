#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 19:50:29 2020

@author: jeffreywang
"""

import argparse
import math 
import random
import sys
import os

parser = argparse.ArgumentParser() 

#parser.add_argument('-rankBy', action = 'store', dest='rank', help= 'rank by: 1 for p adjusted value, 2 for dZsc, 3 for mdist (padj is default)')

parser.add_argument('-geneBed', action = 'store', dest = 'geneBed', help = 'Gene bed for mice/humans')

parser.add_argument('-differentialFile', action = 'store', dest = 'differentialFile', help = 'Differential File: name should be in form: ###_vs_###_full_compartment_details.bedGraph')
                    
parser.add_argument('-cGSEAfile', action = 'store', dest = 'cGSEAfile', help = "Specify GSEA pre-ranked parameters. See documentation for details.")

parser.add_argument('-outputFileName', action = 'store', dest = 'outf', help = "Name of directory for output")

parser.add_argument('-prefix', action = 'store', dest = 'prefix', help = "name of comparison")

results = parser.parse_args()

startdir = os.getcwd()
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))

# Make gene set file
cmd = "Rscript " + os.path.join(scriptdir, "cgsea.r") + " " + results.geneBed + " " + results.differentialFile
os.system(cmd)
geneDiffFile = "Genes." + results.differentialFile

# GSEA specifications
gmtFile = ""
collapse = "No_Collapse"
node = "Max_probe"
norm = "meandiv"
nperm = "1000"
scoring_scheme = "weighted" 
rpt_label = results.prefix
create_svgs = "false"
include_only_symbols = "true"
make_sets = "true"
plot_top_x = "20"
rnd_seed = "timestamp"
set_max = "500"
set_min = "15"
zip_report = "false"
outf = str(results.outf)
ranker = 1

with open(results.cGSEAfile, "r") as params:
    ranker = int(params.readline().strip())
    gmtFile = str(params.readline().strip())
    
#print(type(gmtFile))
class gene:
    def __init__(self, name, value):
        self.name = name
        self.value = value
        self.counter = 1
    def newVal(self, val):
        self.value = min(self.value, val)
        self.counter += 1
    def setVal(self, val):
        self.value = val

genes = []
geneNames = []
vals = []

with open(geneDiffFile, "r") as inf:
    inf.readline()
    for line in inf:
        line = line.strip().split()
        linelength = len(line)
        geneString = line[linelength-1]
        val = 0
        if ranker == 1:
            val = (float(line[linelength-2])) # Padj value = -log10(pval)
        if ranker == 2:
            val = (float(line[linelength-4])) # dZsc
        if ranker == 3:
            val = (float(line[linelength-5])) # M-distance
        geneArr = geneString.split(",")
        #print(geneString)
        for g in geneArr:
            g1 = g.upper() # upper case ~!!
            if g1 not in geneNames:
                geneNames.append(g1)
                newGeneObj = gene(g1, val)
                genes.append(newGeneObj)
            else:
                for a in range(len(genes)):
                    if g1 == genes[a].name:
                        genes[a].newVal(val)

finalGeneArr = genes
finalGeneVals = []
for iterator in range(len(genes)):
    adjval = (-1) * math.log10(finalGeneArr[iterator].value) + 1 
    if adjval in finalGeneVals:
        adjval += random.uniform(0, 0.01)
    else:
        finalGeneVals.append(adjval)
    finalGeneArr[iterator].setVal(adjval)

print(len(finalGeneArr))
finalGeneArr.sort(key=lambda x: x.value, reverse=True)

name = geneDiffFile.split(".")[1].split("_")[0:3]
sep = ""
name = sep.join(name)
print(name + " Top 6 Results")

for a in range(0,6):
    print(finalGeneArr[a].name + " " + str(finalGeneArr[a].value))

name = name + "_genes_ranked.rnk"
with open(name, "w") as outfile:
    for g in finalGeneArr:
        outfile.write(g.name + "\t" + str(g.value) + "\n")
        
if sys.platform.startswith("linux") or sys.platform.startswith("darwin"):
    cmd = os.path.join(scriptdir, "GSEA_4.0.3", "gsea-cli.sh") + " GSEAPreranked -gmx " + gmtFile + " -collapse " + collapse + " -mode " + node + " -norm " + norm + " -nperm " + nperm + " -rnk " + name + " -scoring_scheme " + scoring_scheme + " -rpt_label " + rpt_label + " -create_svgs " + create_svgs + " -include_only_symbols " + include_only_symbols + " -make_sets " + make_sets + " -plot_top_x " + plot_top_x + " -rnd_seed " + rnd_seed + " -set_max " + set_max + " -set_min " + set_min + " -zip_report " + zip_report + " -out " + outf
    print(cmd)
    os.system(cmd)
elif sys.platform.startswith("win32"):
    print(cmd)
    cmd = os.path.join(scriptdir, "GSEA_4.0.3", "gsea-cli.bat") + " GSEAPreranked -gmx " + gmtFile + " -collapse " + collapse + " -mode " + node + " -norm " + norm + " -nperm " + nperm + " -rnk " + name + " -scoring_scheme " + scoring_scheme + " -rpt_label " + rpt_label + " -create_svgs " + create_svgs + " -include_only_symbols " + include_only_symbols + " -make_sets " + make_sets + " -plot_top_x " + plot_top_x + " -rnd_seed " + rnd_seed + " -set_max " + set_max + " -set_min " + set_min + " -zip_report " + zip_report + " -out " + outf
    os.system(cmd)
else:
    print("Current support only for Mac/Windows/Linux.")

# -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /Users/jeffreywang/gsea_home/output/jul03
# gsea-cli.sh GSEAPreranked -gmx /Users/jeffreywang/Downloads/MousePath_All_gmt-Format.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk /Users/jeffreywang/Documents/LJI/ClusterImports/MouseDiff/GSEA_ImportantGenes/NPCvsmESC_genes_ranked.rnk -scoring_scheme weighted -rpt_label GSEA_Mice_ImpGenes -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /Users/jeffreywang/Documents/LJI/ClusterImports/MouseDiff/GSEA_ImportantGenes







