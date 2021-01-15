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

parser.add_argument("-keepIntermediates", action = 'store', dest = 'intermediates', help = "Variable, if set, will keep intermediate files: graphing, R session, tracks, etc.")

parser.add_argument("-blacklist", action = 'store', dest = 'blacklist', help = "Path to bed file with blacklisted regions not to include. Do not include if not wanted.")

parser.add_argument("-ncp", action = 'store', dest = 'ncp', help = "Number of principal components to scan/use. Default is 2.")

parser.add_argument("-repParams", action = 'store', dest = 'repParams', help = "Replicate parameter file (PC variation) to use, if no replicates provided.")

parser.add_argument("-SVfilter", action = 'store', dest = 'filter', help = "SV scores text file. Only for single level HMFA. Chrs should be organized same way as chr.txt.")

parser.add_argument("-removeFile", action = 'store', dest = 'removal', help = "Manual removal of exp + chr. Two columns: Experiment\tChr.")

results = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s : - %(pathname)s - %(levelname)-8s : %(processName)s : %(message)s')

startdir = os.getcwd()
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
svthreshold = 1
if results.goldenpath is not None:
    results.goldenpath = os.path.abspath(results.goldenpath)
if results.blacklist is not None:
    results.blacklist = os.path.abspath(results.blacklist)
if results.repParams is not None:
    results.repParams = os.path.abspath(results.repParams)

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

else:
    print("Missing genome option. System exit.")
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
numGroups = 0
numExp = 0

filteredchrs = {}
humanchrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"] 
micechrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"]
possiblechrs = []
if results.genome.lower() == "mm10" or results.genome.lower() == "mm9":
    possiblechrs = micechrs
elif results.genome.lower() == "hg38" or results.genome.lower() == "hg19":
    possiblechrs = humanchrs
else:
    print("WARNING/ERROR: Not using mice or human dataset.")
    
for elem in possiblechrs:
    chrTag = elem
    filteredchrs[chrTag] = []
    
if results.filter is not None:
    numlines = 0
    with open(results.filter, "r") as filterfile:
        exps_list = filterfile.readline().strip().split()
        for line in filterfile:
            l = line.strip().split()
            chrTag = l[0][l[0].index("chr")+3:]
            badExps = []
            for a in range(1, len(l)):
                if float(l[a]) > 1:
                    badExps.append(exps_list[a-1])
            if chrTag in filteredchrs:     
                for a in badExps:
                    filteredchrs[chrTag].append(a)

if results.removal is not None:
    with open(results.removal, "r") as removefile:
        for line in removefile:
            l = line.strip().split()
            if l[0] not in filteredchrs[l[1]]:
                filteredchrs[l[1]].append(l[0])
    
def getGroups(chrTag):
    global names
    global groups
    global groups_excl
    global group_sizes
    global destinations
    global groupings
    global groupings_excl
    global groupings_sizes
    global isGrouping
    global startdir
    global currGroups
    global filteredChrs
    global numGroups
    global numExp
    
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
            if results.filter is not None:
                if temp[1] in filteredchrs[chrTag]: # could be temp[1]
                    continue
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
                print("WARNING. Using SV filtering WITH two-tier grouping. Not recommended but continuing.")
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
    getGroups(chrNum)
    workdir = startdir
    os.chdir(startdir)
    newdir = "chr_" + chrNum
    if os.path.exists(newdir) is False:
        os.mkdir(newdir)
    os.chdir(newdir)
    chrTag = chrNum    
    cmd = "python " + os.path.join(scriptdir, "run.py") + " -nExp " + str(numExp) + " -chrNum " + chrTag + " -res " + results.res + " -numGroups " + str(numGroups) + " -grouping 1"
    if results.ncp is not None:
        cmd = cmd + " -ncp " + str(results.ncp)
    else:
        cmd = cmd + " -ncp 2"
    for a in group_sizes:
        cmd = cmd + " -group " + str(a)
    for a in names:
        cmd = cmd + " -expNames " + a
    for a in groups_excl:
        cmd = cmd + " -groupNames " + a
    if results.blacklist is not None:
        cmd = cmd + " -blacklist " + results.blacklist
    if isGrouping:
        for g in groupings_excl:
            cmd = cmd + " -groupings " + g
        for g in groupings_sizes:
            cmd = cmd + " -groupingNums " + str(g)
    
    if results.genome is not None:
        cmd = cmd + " -genome " + str(results.genome)
    else:
        print("ERROR: Specify Genome.")
        sys.exit(1)
        
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
                    b = file.find(".") # change made with NoArm in mind
                    chrIden = (file[(a+3):b])
                    if chrIden == chrNum:
                        file_path = file
                        break
            cmd = cmd + " -prePath " + os.path.join(elem, file_path)

    
    print("\n" + cmd + "\n")
    os.system(cmd)    
    os.mkdir("pcFiles")
    incr = 1
    for x in range(1, (len(names)+1)):
        pcaFileLocation = "pc_" + str(names[x-1]) + "_exp_" + str(incr) + ".txt" 
        shutil.move(pcaFileLocation, "pcFiles")
        incr+=1
    
    incr = 1
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
        getGroups(chrNum)
        print(names)
        print(groups)
        newdir = "chr_" + chrNum
        if os.path.exists(newdir) is False:
            os.mkdir(newdir)
        print("\n###################")
        print("Chromosome " + chrNum + " Output")
        print("###################")
        print(scriptdir)
        print(numExp)
        print(numGroups)
        
        chrTag = chrNum
            
        cmd = "python " + os.path.join(scriptdir, "run.py") + " -nExp " + str(numExp) + " -chrNum " + chrTag + " -res " + results.res + " -numGroups " + str(numGroups) + " -grouping 1"
        if results.ncp is not None:
            cmd = cmd + " -ncp " + str(results.ncp)
        else:
            cmd = cmd + " -ncp 2"
        for a in group_sizes:
            cmd = cmd + " -group " + str(a)
        for a in names:
            cmd = cmd + " -expNames " + a
        for a in groups_excl:
            cmd = cmd + " -groupNames " + a
        if results.blacklist is not None:
            cmd = cmd + " -blacklist " + results.blacklist
        if results.genome is not None:
            cmd = cmd + " -genome " + str(results.genome)
        else:
            print("ERROR: Specify Genome.")
            sys.exit(1)
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
        for file in glob.glob("hmfa_*"):
            shutil.move(file, newdir)
        iterator = 1
        for x in range(1, (len(groups_excl)+1)):
            pcaFileLocation = os.path.join(newdir, "hmfa_" + str(groups_excl[x-1]) + "_exp_" + str(x) + ".txt")
            cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -eigfile " + pcaFileLocation + " -chr "+ chrNum + " -exp " + str(groups_excl[x-1])
            print(cmd)
            os.system(cmd)
        if os.path.exists("pcFiles") == False:
            os.mkdir("pcFiles")
        for x in range(1, (len(names)+1)):
            pcaFileLocation = "pc_" + str(names[x-1]) + "_exp_" + str(x) + ".txt" 
            shutil.move(pcaFileLocation, "pcFiles")
        shutil.move("pcFiles", newdir)
        for file in glob.glob("BalancedChrMatrix*"):
            shutil.move(file, newdir)
        for file in glob.glob("*_grouping.txt"):
            shutil.move(file, newdir)
        for file in glob.glob("*.bedGraph"):
            shutil.move(file, newdir)
        for file in glob.glob("coordinates.txt"):
            shutil.move(file, newdir)
        for file in glob.glob("lengthdoc.txt"):
            shutil.move(file, newdir)
        for file in glob.glob("compartmentSwitch_*"):
            shutil.move(file, newdir)
        for file in glob.glob("PCselect*"):
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
    if results.intermediates is None:
        for file in glob.glob("compartmentSwitch_*"):
            os.remove(file)
        for file in glob.glob("Rsession*"):
            os.remove(file)
        for file in glob.glob("*.bed"):
            os.remove(file)
    else:
        os.mkdir("comptSwitch")
        for file in glob.glob("compartmentSwitch_*"):
            shutil.move(file, "comptSwitch")
    os.chdir("..")
    
# =============================================================================
# Exiting Chromosome-By-Chromosome Analysis
# =============================================================================

print("\n###################")
print("FINISHED Chromosome-By-Chromosome Analysis. Now Going to Summary + Visualization.")
print("###################")

# =============================================================================
# PC-Info Document
# =============================================================================

ncp = 2
if results.ncp is not None:
    ncp = int(results.ncp)

with open("chr_info.txt", "w") as info:
    info.write("chr\tpc.select\tsign\tExperiment\n")
    for chrNum in chrlist:
        getGroups(chrNum)
        name = "chr_" + chrNum
        os.chdir(name)
        with open("PCselection.txt") as pcdoc: # pc1.gc  pc2.gc  id      chr     pc.select       sign    pc1.len pc2.len pc.len.warning
            pcdoc.readline()
            it = 0
            numGroups = len(groups)
            for a in range(numGroups):
                for b in range(ncp):
                    pcnum = b+1
                    line = pcdoc.readline()
                    l = line.strip().split()
                    #print(l)
                    if len(l) == 0:
                        continue
                    expname = l[len(l)-2]
                    sign = l[3]
                    if l[len(l)-1] == "yes":
                        info.write(name + "\t" + str(pcnum) + "\t" + sign + "\t" + expname + "\n")
        os.chdir("..")

# # =============================================================================
# # Differential Calling
# # =============================================================================

cmd = "python " + os.path.join(scriptdir, "differentialCalling.py") + " -inputFile " + results.input + " -chrFile " + results.chrs + " -makePlots 1 -res " + str(results.res) + " -genome " + results.genome 
if isGrouping:
    if len(groupings_excl) > 2:
        cmd = cmd + " -multiComp 1"
else:
    if len(groups_excl) > 2:
        cmd = cmd + " -multiComp 1"
if results.filter is not None:
    cmd = cmd + "  -SVfilter " + results.filter 
if results.removal is not None:
    cmd = cmd + " -removeFile " + results.removal 
if results.repParams is not None:
    cmd = cmd + " -repParams " + results.repParams
if results.blacklist is not None:
    cmd = cmd + " -blacklist " + results.blacklist
print(cmd)
os.system(cmd)

# =============================================================================
# Creating score plots
# =============================================================================

if results.filter is None and results.removal is not None:
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
