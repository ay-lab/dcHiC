#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 11:42:17 2020

@author: jeffreywang
"""

import argparse
import os
import shutil
import glob
import multiprocessing
import logging
import sys

parser = argparse.ArgumentParser()

parser.add_argument("-res", action = 'store', dest = 'res', help = "Resolution of files")

parser.add_argument("-inputFile", action = "store", dest = 'input', help = "Input file")

parser.add_argument("-chrFile", action = 'store', dest = 'chrs', help = "file with chromosome specifications: one chromosome per line. Must be entered in order.")

parser.add_argument("-input", action = 'store', dest = 'inputNum', help = "Input format (for naming): HOMER is 1, all else is 2")

parser.add_argument("-parallel", action = 'store', dest = 'par', help = "Use if parallel processing is wanted (otherwise will do chromosomes in sequence). Enter number of cores to be used or enter 0 to use max number of cores possoble. Using parallel implementation will cause the output to be messy.")

parser.add_argument("-alignData", action = 'store', dest = 'goldenpath', help = "[REQUIRED FOR PARALLEL PROCESSING] Absolute location of GoldenPath data")

parser.add_argument('-genome', action = 'store', dest = 'genome', help = "Used to determine eigenvector sign. Options: hg38, hg19, mm10, mm9 [If not provided, analysis will skip this step]")

parser.add_argument('-signAnalysis', action = 'store', dest = 'analysis', help = "Used to determine which data will be used for determining eigenvector sign. Options: \"tss\", \"gc\". [If not provided, analysis will skip this step.]")

parser.add_argument('-cGSEA', action = 'store', dest = 'cgsea', help = "Set if you want GSEA performed. File with basic cGSEA specifications as specified in documentation... Do not include field if cGSEA not wanted.")

parser.add_argument("-keepIntermediates", action = 'store', dest = 'intermediates', help = "Variable, if set, will keep intermediate files: graphing, R session, tracks, etc.")

results = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s : - %(pathname)s - %(levelname)-8s : %(processName)s : %(message)s')

startdir = os.getcwd()
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
    
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

tempchrlist = []
with open(results.chrs, "r") as input:
    for line in input:
        a = line.strip()
        tempchrlist.append(a)

chrlist = [] # Remove Duplicates, if any
for chrelem in tempchrlist:
    if chrelem not in chrlist:
        chrlist.append(chrelem)
print(chrlist)
     
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

# =============================================================================
# Method for parallel processing; called later
# =============================================================================

def chr_process(chrNum):
    workdir = startdir
    os.chdir(startdir)
    newdir = "chr_" + chrNum
    cmd = "mkdir " + newdir
    os.system(cmd)
    os.chdir(newdir)
    cmd = "python " + os.path.join(scriptdir, "run.py") + " -nExp " + str(numExp) + " -chrNum " + chrNum + " -res " + results.res + " -numGroups " + str(numGroups) + " -grouping 1"
    for a in group_sizes:
        cmd = cmd + " -group " + str(a)
    for a in names:
        cmd = cmd + " -expNames " + a
    for a in groups_excl:
        cmd = cmd + " -groupNames " + a
    if isGrouping:
        for g in groupings_excl:
            cmd = cmd + " -groupings " + g
        for g in groupings_sizes:
            cmd = cmd + " -groupingNums " + str(g)
    
    if results.genome is not None and results.analysis is not None:
        cmd = cmd + " -genome " + str(results.genome) + " -signAnalysis " + str(results.analysis)
    
    if results.goldenpath is not None:
        cmd = cmd + " -alignData " + str(results.goldenpath)
        
    if int(results.inputNum) == 1: # HOMER generated matrices
        for elem in destinations:
            file_list = os.listdir(elem)
            file_path = ""
            for file in file_list:
                a = file.find("chr", 0)
                b = file.find(".", 3)
                chrIden = (file[(a+3):b])
                if chrIden == chrNum:
                    file_path = file
                    break
            cmd = cmd + " -prePath " + os.path.join(elem, file_path)
            
    if int(results.inputNum) == 2: # program-generated matrices
        for elem in destinations:
            file_list = os.listdir(elem)
            file_path = ""
            for file in file_list:
                if "matrix" in file:
                    a = file.find("chr")
                    b = file.find(".")
                    chrIden = (file[(a+3):b])
                    if chrIden == chrNum:
                        file_path = file
                        break
            cmd = cmd + " -prePath " + os.path.join(elem, file_path)
    
    print("\n" + cmd + "\n")
    os.system(cmd)    
    os.mkdir("pcFiles")
    for x in range(1, (len(names)+1)):
        pcaFileLocation = "pc1_" + str(names[x-1]) + "_exp_" + str(x) + ".txt" 
        shutil.move(pcaFileLocation, "pcFiles")
        pcaFileLocation = "pc2_" + str(names[x-1]) + "_exp_" + str(x) + ".txt" 
        shutil.move(pcaFileLocation, "pcFiles")
        
    for x in range(1, (len(groups_excl)+1)):
        pcaFileLocation = "hmfa_" + str(groups_excl[x-1]) + "_exp_" + str(x) + ".txt" # this will have to change if naming conventions change
        cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -eigfile " + pcaFileLocation + " -chr "+ chrNum + " -exp " +  str(groups_excl[x-1])
        logging.debug(cmd)
        os.system(cmd)  

    if results.intermediates is None: # this is run in the chromosome directory
        for file in glob.glob("Rsession*"):
            os.remove(file)  
            
    os.chdir(workdir)

# =============================================================================
# Non-parallel processing (chromosomes in sequence)
# =============================================================================

if results.par is None:
    for chrNum in chrlist:
        newdir = "chr_" + chrNum
        if os.path.exists(newdir) is False:
            os.mkdir(newdir)
        print("\n###################")
        print("Chromosome " + chrNum + " Output")
        print("###################")
        print(scriptdir)
        print(numExp)
        print(numGroups)
        cmd = "python " + os.path.join(scriptdir, "run.py") + " -nExp " + str(numExp) + " -chrNum " + chrNum + " -res " + results.res + " -numGroups " + str(numGroups) + " -grouping 1"
        for a in group_sizes:
            cmd = cmd + " -group " + str(a)
        for a in names:
            cmd = cmd + " -expNames " + a
        for a in groups_excl:
            cmd = cmd + " -groupNames " + a
        if results.genome is not None and results.analysis is not None:
            cmd = cmd + " -genome " + str(results.genome) + " -signAnalysis " + str(results.analysis)
        if results.goldenpath is not None:
            cmd = cmd + " -alignData " + str(results.goldenpath)
        if isGrouping:
            for g in groupings_excl:
                cmd = cmd + " -groupings " + g
            for g in groupings_sizes:
                cmd = cmd + " -groupingNums " + str(g)
        if int(results.inputNum) == 1:
            for elem in destinations:
                file_list = os.listdir(elem)
                file_path = ""
                for file in file_list:
                    a = file.find("chr", 0)
                    b = file.find(".", 3)
                    chrIden = (file[(a+3):b])
                    if chrIden == chrNum:
                        file_path = file
                        break
                cmd = cmd + " -prePath " + os.path.join(elem, file_path)
            
        if int(results.inputNum) == 2:
            for elem in destinations:
                file_list = os.listdir(elem)
                file_path = ""
                for file in file_list:
                    if "matrix" in file:
                        a = file.find("chr")
                        b = file.find(".")
                        chrIden = (file[(a+3):b])
                        if chrIden == chrNum:
                            file_path = file
                            break
                cmd = cmd + " -prePath " + os.path.join(elem, file_path)
        print("\n" + cmd + "\n")
        os.system(cmd)
        #for file in glob.glob("BalancedChr*"):
        #    exp_num = file.split("_")[2].split(".")[0]
        #    newname = file.split("_")[0] + "_exp_" + names[int(exp_num)-1] + ".txt"
        #    command = "mv " + file + " " + newname
        #    os.system(command)
        #    shutil.move(newname, newdir)
        for file in glob.glob("hmfa_*"):
            shutil.move(file, newdir)
        iterator = 1
        for x in range(1, (len(groups_excl)+1)):
            pcaFileLocation = os.path.join(newdir, "hmfa_" + str(groups_excl[x-1]) + "_exp_" + str(x) + ".txt")
            cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -eigfile " + pcaFileLocation + " -chr "+ chrNum + " -exp " + str(x)
            print(cmd)
            os.system(cmd)
        if os.path.exists("pcFiles") == False:
            os.mkdir("pcFiles")
        for x in range(1, (len(names)+1)):
            pcaFileLocation = "pc1_" + str(names[x-1]) + "_exp_" + str(x) + ".txt" 
            shutil.move(pcaFileLocation, "pcFiles")
            pcaFileLocation = "pc2_" + str(names[x-1]) + "_exp_" + str(x) + ".txt" 
            shutil.move(pcaFileLocation, "pcFiles")
        shutil.move("pcFiles", newdir)
        for file in glob.glob("*_grouping.txt"):
            shutil.move(file, newdir)
        for file in glob.glob("*.pdf"):
            shutil.move(file, newdir)
        for file in glob.glob("*.bedGraph"):
            shutil.move(file, newdir)
        for file in glob.glob("coordinates.txt"):
            shutil.move(file, newdir)
        for file in glob.glob("lengthdoc.txt"):
            shutil.move(file, newdir)
        for file in glob.glob("compartmentSwitch_*"):
            shutil.move(file, newdir)
        for file in glob.glob("pc_decision*"):
            shutil.move(file, newdir)
        if results.intermediates is None:
            for file in glob.glob("Rsession*"):
                os.remove(file)
        else:
            for file in glob.glob("Rsession*"):
                shutil.move(file, newdir)

# =============================================================================
# Parallel Processing
# =============================================================================
            
else:
    tmp = int(results.par)
    numCores = multiprocessing.cpu_count()
    if tmp > 0:
        numCores = tmp
    logging.debug("\nThere are %d CPUs on this machine" % multiprocessing.cpu_count() + "\n")
    total_tasks = len(chrlist)
    pool = multiprocessing.Pool(numCores)
    res_parallel = pool.map_async(chr_process, chrlist)
    pool.close()
    pool.join()
    logging.debug("Pooling complete.")

# =============================================================================
# Keep/Remove Intermediates 
# =============================================================================

for chrNum in chrlist:
    print(os.getcwd())
    os.chdir("chr_" + chrNum)
    # for file in glob.glob("HMFA*"): # this goes one at a time
    #     exp_num = file.split("_")[3].split(".")[0]
    #     newname = file.split("_")[0] + "_chr" + chrNum + "_" + groups_excl[int(exp_num)-1] + ".bedGraph"
    #     command = "mv " + file + " " + newname
    #     os.system(command)
    if results.intermediates is None:
        for file in glob.glob("compartmentSwitch_*"):
            os.remove(file)
        for file in glob.glob("Rsession*"):
            os.remove(file)
        for file in glob.glob("*.bed"):
            os.remove(file)
    os.chdir("..")
    
# =============================================================================
# Exiting Chromosome-By-Chromosome Analysis
# =============================================================================

print("\n###################")
print("FINISHED Chromosome-By-Chromosome Analysis. Now Going to Summary + Visualization.")
print("###################")

# =============================================================================
# Differential Calling
# =============================================================================

cmd = "python " + os.path.join(scriptdir, "differentialCalling.py") + " -inputFile " + results.input + " -chrFile " + results.chrs + " -makePlots 1 -res " + str(results.res)
if isGrouping:
    if len(groupings_excl) > 2:
        cmd = cmd + " -multiComp 1"
else:
    if len(groups_excl) > 2:
        cmd = cmd + " -multiComp 1"
print(cmd)
os.system(cmd)

# =============================================================================
# Compartmental Gene Set Enrichment Analysis
# =============================================================================

if results.cgsea is not None:
    os.chdir("DifferentialCompartment")
    for diffFile in glob.glob("*full_compartment_details*"):
        if "Genes." in diffFile:
            continue
        geneBedName = results.genome.lower() + "_gene_pos.bed"
        pos = diffFile.index("full_compartment_details") -1
        name = diffFile[:pos]
        cmd = "python " + os.path.join(scriptdir, "cgsea.py") + " -differentialFile " + diffFile + " -cGSEAfile ../" + results.cgsea + " -prefix " + name + " -outputFileName GSEA_" + name
        print(cmd)
        os.system(cmd)
    os.chdir("..")
        
# =============================================================================
# Creating score plots
# =============================================================================

command = "python " + os.path.join(scriptdir, "coordinates.py") + " -nExp " + str(numExp) + " " # coordinates averaging
for a in names:
    command += "-exp " + a + " "
command += "-nGroups " + str(len(groups_excl)) + " "
for x in groups_excl: # should be groups_excl
    command += "-groups " + x + " "
for x in chrlist:
    command += "-chr " + x + " "
print("\n" + command + "\n")
os.system(command)

# =============================================================================
# PC-Info Document
# =============================================================================

with open("chr_info.txt", "w") as info:
    for chrNum in chrlist:
        os.chdir("chr_" + chrNum)
        pcNum = 0
        with open("lengthdoc.txt", "r") as lengthdoc:
            a = lengthdoc.readline()
            b = lengthdoc.readline()
            pcNum = lengthdoc.readline()
        info.write("Chromosome " + str(chrNum) + ": PC " + str(pcNum)+"\n")
        os.chdir("..")
