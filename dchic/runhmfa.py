#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 11:25:36 2020

@author: jeffreywang
"""

import argparse
import os
import glob
import shutil
import sys

parser = argparse.ArgumentParser()

parser.add_argument('-res', action = 'store', dest = 'res', help = 'resolution (10000, 25000, 50000, etc.)')

parser.add_argument("-inputFile", action = "store", dest = 'input', help = "Input file")

parser.add_argument('-chr', action = 'store', dest = 'chr', help = 'chromosome number: 1, 2, 3...X, Y')

parser.add_argument('-genome', action = 'store', dest = 'genome', help = "Used to determine eigenvector sign. Options: hg38, hg19, mm10, mm9 [If not provided, analysis will skip this step]")

parser.add_argument('-signAnalysis', action = 'store', dest = 'analysis', help = "Used to determine which data will be used for determining eigenvector sign. Options: TSS, GC content. [If not provided, analysis will skip this step.]")

parser.add_argument("-alignData", action = 'store', dest = 'goldenpath', help = "Absolute location of GoldenPath data [If not provided, analysis will skip this step.]")

results = parser.parse_args()

startdir = os.getcwd()
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
chrNum = str(results.chr)
    
# =============================================================================
# System Checks
# =============================================================================

if results.genome is not None: 
    validGenomes = ["hg38", "hg19", "mm10", "mm9"]
    valid = False
    for v in validGenomes:
        if results.genome == v:
            valid = True
            break
    if valid == False:
        print("Bad genome entered. Try again. System exit.")
        sys.exit(1)

if results.analysis is not None:
    validAnalyses = ["tss", "gc"]
    valid = False
    for v in validAnalyses:
        if results.analysis == v:
            valid = True
            break
    if valid == False:
        print("Bad analysis method entered. Try again. System exit.")
        sys.exit(1)

# =============================================================================
# Scanning input
# =============================================================================
        
names = []
groups = [] # layer 1 of organization 
groups_excl = []
group_sizes = []
destinations = []

groupings = [] # groupings of groups
groupings_excl = []
groupings_sizes = []
isGrouping = False
startdir = os.getcwd()
currGroups = []

with open(results.input, 'r') as input:
    for line in input:
        if line == "\n":
            continue
        line = line.strip()
        temp = line.split()
        names.append(temp[0])
        groups.append(temp[1])
        if temp[1] not in groups_excl:
            groups_excl.append(temp[1])
            group_sizes.append(0)
        group_sizes[groups_excl.index(temp[1])] += 1
        if len(temp) == 3:
            destinations.append(temp[2])
        elif len(temp) == 4:
            isGrouping = True
            destinations.append(temp[3])
            groupings.append(temp[2])
            if temp[2] not in groupings_excl:
                currGroups.clear()
                currGroups.append(temp[1])
                groupings_excl.append(temp[2])
                groupings_sizes.append(0)
                groupings_sizes[groupings_excl.index(temp[2])] += 1
            if temp[1] not in currGroups:
                groupings_sizes[groupings_excl.index(temp[2])] += 1
        else:
            print("Error in input file formatting. Exiting.")
            sys.exit(1)

numGroups = len(groups_excl)
numExp = len(names)

cmd = "python " + os.path.join(scriptdir, "run.py") + " -nExp " + str(len(destinations)) + " -chrNum " + chrNum + " -res " + results.res + " -numGroups " + str(numGroups) + " -grouping 1" # change to 1 or 2
for a in group_sizes:
    cmd = cmd + " -group " + str(a)
for a in names:
    cmd = cmd + " -expNames " + a
for a in groups_excl:
    cmd = cmd + " -groupNames " + a
for a in destinations:
    cmd = cmd + " -prePath " + a
if results.genome is not None and results.analysis is not None:
    cmd = cmd + " -genome " + str(results.genome) + " -signAnalysis " + str(results.analysis)
if results.goldenpath is not None:
    cmd = cmd + " -alignData " + str(results.goldenpath)
if isGrouping:
    for g in groupings_excl:
        cmd = cmd + " -groupings " + g
    for g in groupings_sizes:
        cmd = cmd + " -groupingNums " + str(g)

print("\n" + cmd + "\n")
os.system(cmd)

for x in range(1, (len(groups_excl)+1)):
    pcaFileLocation = "hmfa_" + str(groups_excl[x-1]) + "_exp_" + str(x) + ".txt" # this will have to change if naming conventions change
    cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -eigfile " + pcaFileLocation + " -chr "+ chrNum + " -exp " + str(groups_excl[x-1])
    print(cmd)
    os.system(cmd)

#for file in glob.glob("HMFA*"): # this goes one at a time
#    exp_num = file.split("_")[3].split(".")[0]
#    newname = file.split("_")[0] + "_chr" + chrNum + "_" + groups_excl[int(exp_num)-1] + ".bedGraph"
#    command = "mv " + file + " " + newname
#    os.system(command)

os.mkdir("pcFiles")
for x in range(1, (len(names)+1)):
    pcaFileLocation = "pc1_" + str(names[x-1]) + "_exp_" + str(x) + ".txt" 
    shutil.move(pcaFileLocation, "pcFiles")
    pcaFileLocation = "pc2_" + str(names[x-1]) + "_exp_" + str(x) + ".txt" 
    shutil.move(pcaFileLocation, "pcFiles")

os.mkdir("comptSwitch")
for a in glob.glob("compartmentSwitch*"):
    shutil.move(a, "comptSwitch")
