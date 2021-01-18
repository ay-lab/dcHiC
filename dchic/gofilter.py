#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:12:14 2020

@author: jeffreywang
"""

import argparse
import os
import sys
import ast

parser = argparse.ArgumentParser()
parser.add_argument("-dir", action = 'store', dest = 'dir', help = "Abs. path of overhead directory of all results")
parser.add_argument("-diffcompt", action = 'store', dest = 'diffcompt', help = "Differential Compartment File")
parser.add_argument("-config", action = 'store', dest = 'config', help = "Configuration file. Line 1 = include; Line 2 = exclude. Comma separate exp. names with no spaces.")
parser.add_argument("-outprefix", action = 'store', dest = 'outputdir', help = "Name of output directory and final file prefix")
parser.add_argument("-genome", action = 'store', dest = 'genome', help = "genome: mm10/mm9/hg38/hg19")
parser.add_argument("-geneBed", action = 'store', dest = 'geneBed', help = "gene bed for genome")
parser.add_argument("-runOption", action = 'store', dest = 'runOp', help = "1 for direction, 2 for sign, 3 for both")
parser.add_argument("-orientation", action = 'store', dest = 'orient', help = "1 = Active to Inactive; 2 = Inactive to Active. Default for enrichment is 1.")
parser.add_argument("-slack", action = 'store', dest = 'slack', help = "# of bins (surrounding differential region) to count. Default = 0.")
results = parser.parse_args()

scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
results.config = os.path.abspath(results.config)
results.diffcompt = os.path.abspath(results.diffcompt)
results.geneBed = os.path.abspath(results.geneBed)
results.dir = os.path.abspath(results.dir)
if os.path.isdir(results.outputdir) == False:
    os.mkdir(results.outputdir)
os.chdir(results.outputdir)
startdir = os.getcwd()
slack = 0
if results.slack is not None:
    slack = int(results.slack)

print("Slack Given In GO Analysis: " + str(slack))
    
resultsdir = os.path.abspath(results.dir)

orientation = "Active2Inactive"
if results.orient is not None:
    if results.orient == "1":
        orientation = "Active2Inactive"
    elif results.orient == "2":
        orientation = "Inactive2Active"
else:
    orientation = "Active2Inactive"

if results.runOp == "1":
    runOp = "direction"
elif results.runOp == "2":
    runOp = "sign"
elif results.runOp == "3":
    runOp = "both"
else:
    print("Sign option mode entered incorrectly. Chose 1, 2, or 3. Exit")
    sys.exit(1)
    
with open(results.config, "r") as config:
    include = config.readline().strip()
    exclude = config.readline().strip()

cmd = "Rscript " + os.path.join(scriptdir, "cluster.r") + " " + runOp + " " + orientation + " " + exclude + " " + include + " " + resultsdir + " 1" 
print(cmd)
os.system(cmd)

diff = results.diffcompt
tmpfile = "tmp"
finalfile = results.outputdir + ".bedGraph"

diffbins = []
diffbinsitems = []
init = ""
with open(diff, "r") as inf:
    init = inf.readline()
    for line in inf:
        l = line.strip().split()
        tag = l[0] + "\t" + l[1] + "\t" + l[2]
        diffbins.append(tag)
        diffbinsitems.append(l[3:])

tmpbins = []
tmpbinitems = []
with open(tmpfile, "r") as inf:
    inf.readline()
    for line in inf:
        l = line.strip().split()
        tag = l[0] + "\t" + l[1] + "\t" + l[2]
        tmpbins.append(tag)
        tmpbinitems.append(l[3:])

intersection = []

for a in range(len(diffbins)):
    if diffbins[a] in tmpbins:
        intersection.append(a)

with open(finalfile, "w") as file:
    file.write(init)
    for b in range(len(diffbinsitems)):
        if b in intersection:
            file.write(diffbins[b] + "\t" + "\t".join(diffbinsitems[b]) + "\n")

cmd = "Rscript " + os.path.join(scriptdir, "cluster.r") + " " + str(slack) + " " + results.geneBed + " " + finalfile + " 2"
print(cmd)
os.system(cmd)

genelist = []
with open("Genes." + finalfile, "r") as inf:
    for line in inf:
        l = line.strip().split()
        geneArr = l[len(l)-1].split(",")
        for a in geneArr:
            gene = a.lower()
            if gene not in genelist:
                genelist.append(gene)

if len(genelist) == 0:
    print("No genes found. Stopping")
    sys.exit(1)

else:
    entrezlist = []
    finalfile = "GO_" + results.outputdir + ".txt"
    cmd = "curl -H 'Content-Type: text/json' -d '{\"Symbols\":[\"" + genelist[0] + "\""
    for a in range(1, len(genelist)):
        cmd = cmd + ",\"" + genelist[a] + "\""
    cmd = cmd + "]}' https://toppgene.cchmc.org/API/lookup"  
    print("\n"+cmd+"\n")
    res = os.popen(cmd).read()
    res = res[res.index("["):].strip()
    res = res[:len(res)-1]
    results = ast.literal_eval(res)
    for a in results:
        entrezlist.append(a['Entrez'])
    passin = "{\"Genes\":[" + str(entrezlist[0])
    for a in range(1, len(entrezlist)):
        passin = passin + "," + str(entrezlist[a])
    passin = passin + "],\"Categories\":[{\"Type\": \"GeneOntologyBiologicalProcess\", \"PValue\": 0.05, \"MinGenes\": 1, \"MaxGenes\": 1500, \"MaxResults\": 30, \"Correction\": \"FDR\"}]}"
    cmd = "curl -H 'Content-Type: text/json' -d '" + passin + "' https://toppgene.cchmc.org/API/enrich"
    print("\n"+cmd+"\n")
    res = os.popen(cmd).read()
    if "null" in res and len(res) < 25:
        print("No Significant GO Annotations Found. Exiting.")
        sys.exit(0)
    res = res[:len(res)-1]
    results = ast.literal_eval(res)
    with open(finalfile, "w") as outf:
        outf.write("GO ID\tGO Term\tP-value\tFDR B-H\tFDR B-Y\tFDR Bonferonni\tTotal Genes\tGenes In Term\tGenes In Query\tGenes In Term & Query\tGenes\n")
        for a in range(len(results)):
            GeneList = results[a]['Genes'][0]['Symbol']
            for b in range(1, len(results[a]['Genes'])):
                GeneList = GeneList + "," + results[a]['Genes'][b]['Symbol']
            outstr = str(results[a]['ID']) + "\t" + str(results[a]['Name']) + "\t" + str(results[a]['PValue']) + "\t" + str(results[a]['QValueFDRBH']) + "\t" + str(results[a]['QValueFDRBY']) + "\t" + str(results[a]['QValueBonferroni']) + "\t" + str(results[a]['TotalGenes']) + "\t" + str(results[a]['GenesInTerm']) + "\t" + str(results[a]['GenesInTermInQuery']) + "\t" + str(GeneList) + "\n" 
            outf.write(outstr)
    
