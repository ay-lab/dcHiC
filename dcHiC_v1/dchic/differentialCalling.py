#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 18:28:11 2020

@author: jeffreywang
"""
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
import math
import statistics
from statistics import stdev
from scipy import stats
import shutil
import sys

parser = argparse.ArgumentParser()

parser.add_argument("-inputFile", action = 'store', dest = 'inputfile', help = "input text file (input.txt)")

parser.add_argument("-chrFile", action = 'store', dest = 'chr', help = "Chr file, same as one used in dchic.py")

parser.add_argument("-makePlots", action = 'store', dest = 'makePlots', help = "Set to true if you want replicate variation plots to be made (basically visualizes how much replicate variation there is). DEBUG feature. Set true by default.")

parser.add_argument("-res", action = 'store', dest = 'res', help = "Resolution")

parser.add_argument("-multiComp", action = 'store', dest = 'mcomp', help = "Include this field (-multiComp 1) if you want a combined comparison. DEBUG feature. Set true by default.")

parser.add_argument("-repParams", action = 'store', dest = 'reps', help = "Replicate parameter file (PC variation) to use, if no replicates provided.")

parser.add_argument("-blacklist", action = 'store', dest = 'blacklist', help = "Blacklist, if used.")

parser.add_argument("-genome", action = 'store', dest = 'genome', help = "Genome: hg38/19, mm10/9")

#parser.add_argument("-SVfilter", action = 'store', dest = 'filter', help = "DEBUG feature in development. SV scores text file. Only for single level HMFA. Chrs should be organized same way as chr.txt.")

#parser.add_argument("-removeFile", action = 'store', dest = 'removal', help = "DEBUG feature in development. Manual removal of exp + chr. Two columns: Experiment\tChr.")

#parser.add_argument("-keepIntermediates", action = 'store', dest = "keepIntermediates", help = "DEBUG FEATURE: Activate to output replicate fitting info while processing")

results = parser.parse_args()

startdir = os.getcwd()
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))

if results.blacklist is not None:
    results.blacklist = os.path.abspath(results.blacklist)

pc_coords = []
bins = []

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

filterStatus = 0 # results.filter == none --> 0
removeStatus = 0 # results.removal == none --> 0

tempchrlist = []
with open(results.chr, "r") as input:
    for line in input:
        a = line.strip()
        tempchrlist.append(a)

chrlist = [] # Remove Duplicates, if any
for chrelem in tempchrlist:
    if chrelem not in chrlist:
        chrlist.append(chrelem)

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
    
if filterStatus != 0:
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

if removeStatus != 0:
    with open(results.removal, "r") as removefile:
        for line in removefile:
            l = line.strip().split()
            if l[0] not in filteredchrs[l[1]]:
                filteredchrs[l[1]].append(l[0])

if isGrouping == False:
    if (len(names) != len(groups)):
        print("Erroneous input file. Length of names, groups don't match.")
else:
    if (len(names) != len(groups)) or (len(groups) != len(groupings)):
        print("Erroneous input file. Length of names, groups, and groupings don't match.")

def needsSeparateDiffCmp():
    global chrlist
    totalchrs = len(chrlist)
    affected = 0
    for c in chrlist:
        if len(filteredchrs[c]) > 0:
            affected += 1
    if float(affected/totalchrs) > 0.2:
        return True
    else:
        return False

def getUnaffectedChrs():
    global chrlist
    goodlist = []
    badlist = []
    for c in chrlist:
        if len(filteredchrs[c]) == 0:
            goodlist.append(c)
        else:
            badlist.append(c)
    re = [goodlist, badlist]
    return re
        
def setDestinations(Chr):
    global destinations
    destinations = []
    currdir = os.getcwd()
    chrTag = "chr_" + Chr
    pcFiles = os.path.join(currdir, chrTag, "pcFiles")
    for a in range(len(names)):
        filename = "pc_" + str(names[a]) + "_exp_" + str(a+1) + ".txt"
        destinations.append(os.path.join(pcFiles, filename))
    
    # pc_decision_file = os.path.join(currdir, chrTag, "pc_decision.txt")
    # with open(pc_decision_file, "r") as pc_file:
    #     pc_decision = int(pc_file.readline().strip())
    # if pc_decision == 1:
    #     for a in range(len(names)):
    #         filename = "pc1_" + str(names[a]) + "_exp_" + str(a+1) + ".txt"
    #         destinations.append(os.path.join(pcFiles, filename))
    # if pc_decision == 2:
    #     for a in range(len(names)):
    #         filename = "pc2_" + str(names[a]) + "_exp_" + str(a+1) + ".txt"
    #         destinations.append(os.path.join(pcFiles, filename))

def checkInputs(chrTag):
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
    
    with open(results.inputfile, 'r') as input:
        for line in input:
            line = line.strip()
            temp = line.split()
            if filterStatus != 0 and chrTag != 0:
                if temp[1] in filteredchrs[chrTag]:
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
                destinations.append(temp[3])
                groupings.append(temp[2])
                if temp[2] not in groupings_excl:
                    groupings_excl.append(temp[2])
                    groupings_sizes.append(0)
                groupings_sizes[groupings_excl.index(temp[2])] += 1
            else:
                print("Error in input file formatting. Exiting.")
                sys.exit(1)

def makeSampleFile():
    global names
    global groups
    global groupings
    print("Sample File in Progress")
    #print(names)
    #print(groups)
    for a in range(len(names)):
        if "-" in names[a]:
            names[a] = names[a].replace("-", ".")
    for a in range(len(groupings)):
        if "-" in groupings[a]:
            groupings[a] = groupings[a].replace("-", ".")
    for a in range(len(groups)):
        if "-" in groups[a]:
            groups[a] = groups[a].replace("-", ".")
            
    with open("samplefile.txt", "w") as outf:
        if isGrouping:
            outf.write("replicate\tprefix\tgroup\n")
            for a in range(len(names)):
                outf.write(names[a] +"\t" + groups[a] + "\t" + groupings[a] + "\n")
        else:
            outf.write("replicate\tprefix\n")
            for a in range(len(names)):
                outf.write(names[a] +"\t" + groups[a] + "\n")
                
def getData(file1, file2):
    r1data = []
    r2data = []
    with open(file1, "r") as r1:
        r1.readline()
        for line in r1:
            linearr = line.split()
            r1data.append(float(linearr[1]))
    
    with open(file2, "r") as r2:
        r2.readline()
        for line in r2:
            linearr = line.split()
            r2data.append(float(linearr[1]))
            bins.append(linearr[0])
    
    box = [r1data, r2data]
    return box

def createModel(r1data, r2data, chrNum, makeImages):
    x = np.array(r1data).reshape((-1, 1))
    y = np.array(r2data)
    model = LinearRegression().fit(x, y)
    if makeImages: 
        plt.plot(x, y, 'ro', markersize = 1)
        plt.plot(x, model.predict(x))
        plt.xlabel("Replicate PC Values")
        plt.ylabel("Replicate PC Values")
        #plt.show()
        filename = "Chr" + str(chrNum) + "_replicates.png"
        plt.savefig(filename, bbox_inches='tight', dpi = 600)
        plt.close()
    mod1 = model.coef_[0]
    mod2 = model.intercept_
    arr = [mod1, mod2] # [slope, intercept]
    return arr

def pointLineDist(p1x, p1y, slope, intercept):
    val = abs((p1x * slope * -1) + (p1y) - intercept) / math.sqrt(slope ** 2 + 1)
    return val

def makeDirectory():
    if os.path.isdir("ReplicateImages") == False:
        os.mkdir("ReplicateImages")

def getRepParams(chrNum): 
    global chrlist
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
    
    if results.makePlots is not None:
        makeDirectory()
    if isinstance(chrNum, list):
        chrlist = chrNum
    elif chrNum != 0: # remove / filter option in place
        checkInputs(chrNum)
        chrlist = [chrNum]
    
    #print(destinations)
    with open("chrdistances.txt", "w") as eucDistanceFile:
        eucDistanceFile.write("m\ts\tchr\n")
        for Chr in chrlist:
            r1data = []
            r2data = []
            setDestinations(Chr)
            #shutil.copy("pcfiles.txt", "chr_" + Chr)
            distances = []
            chrTag = "Chr" + str(Chr)
            for group in range(len(groups_excl)):  
                size = group_sizes[group]
                if size < 2: # Change made to accomodate some non-replicate data
                    continue
                destpos = 0
                for c in range(group):
                    destpos += group_sizes[c]
                for a in range(size):
                    for b in range(a+1, size):
                        box = getData(destinations[a + destpos], destinations[b+destpos])
                        r1data.extend(stats.zscore(box[0], ddof = 1).tolist())
                        r2data.extend(stats.zscore(box[1], ddof = 1).tolist())
            if results.makePlots is not None:
                vals = createModel(r1data, r2data, str(Chr), True) 
                pngName = chrTag + "_replicates.png"
                shutil.move(pngName, "ReplicateImages")
            else:
                vals = createModel(r1data, r2data, str(Chr), False)
            for inc in range(len(r1data)):
                dist = (pointLineDist(r1data[inc], r2data[inc], vals[0], vals[1]))
                distances.append(dist)
            chrMean = statistics.mean(distances)
            chrStd = stdev(distances)
            print("Chromosome " + Chr + " /// Mean: " + str(chrMean) + " Std: " + str(chrStd))
            eucDistanceFile.write(str(chrMean) + "\t" + str(chrStd) + "\t" + str(Chr) + "\n")
    

def main():
    if results.reps is not None:
        distfile = results.reps
    else:
        distfile = "chrdistances.txt"
    if filterStatus == 0 and removeStatus == 0: # run all combined
        checkInputs(0)
        if results.reps is None:
            getRepParams(0)
        makeSampleFile()
        for Chr in chrlist:
            chrtag = "chr_" + Chr
            cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -chr " + Chr + " -combined 1 -res " + results.res
            if results.blacklist is not None:
                cmd = cmd + " -blacklist " + results.blacklist
            print(cmd)
            os.system(cmd)
            
        cmd = "Rscript " + os.path.join(scriptdir, "diffcmp_pythonV.r")
        for Chr in chrlist:
            chrtag = "chr_" + Chr
            cmd = cmd + " " + chrtag + "/" + "chr" + Chr + ".PC.coordinates.txt"
        if results.mcomp is not None:
            cmd = cmd + " 2 "
        else:
            cmd = cmd + " 1 "
        cmd = cmd + " " + results.chr + " " + str(results.res) + " " + distfile + " samplefile.txt " + scriptdir
        print(cmd)
        os.system(cmd)
    
    else:
        if needsSeparateDiffCmp(): # separate DiffCmp per chromosome
            print("NOTE: Because more than 1/5 of chromosomes supplied have a removal, a separate DiffCmp will be run for each one.")
            for Chr in chrlist:
                checkInputs(Chr)
                if results.reps is None:
                    getRepParams(Chr)
                makeSampleFile()
                chrtag = "chr_" + Chr
                cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -chr " + Chr + " -combined 1 -res " + results.res
                if results.blacklist is not None:
                    cmd = cmd + " -blacklist " + results.blacklist
                print(cmd)
                os.system(cmd)
                
                #os.mkdir("chr" + Chr + "_diffcomp")
                cmd = "Rscript " + os.path.join(scriptdir, "diffcmp_pythonV.r") 
                cmd = cmd + " " + chrtag + "/" + "chr" + Chr + ".PC.coordinates.txt"
                if results.mcomp is not None:
                    cmd = cmd + " 2 "
                else:
                    cmd = cmd + " 1 " 
                cmd = cmd + " " + Chr + " " + str(results.res) + " " + distfile + " samplefile.txt " + scriptdir
                print(cmd)
                os.system(cmd)  
        else: # will run one large DiffCmp with non-affected and separate chrwise DiffCmp's for chromosomes with removed cell lines
            relist = getUnaffectedChrs()
            goodlist = relist[0]
            checkInputs(0)
            if results.reps is None:
                getRepParams(goodlist)
            makeSampleFile()
            for Chr in goodlist:
                chrtag = "chr_" + Chr
                cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -chr " + Chr + " -combined 1 -res " + results.res
                if results.blacklist is not None:
                    cmd = cmd + " -blacklist " + results.blacklist
                print(cmd)
                os.system(cmd)
            
            cmd = "Rscript " + os.path.join(scriptdir, "diffcmp_pythonV.r")
            for Chr in goodlist:
                chrtag = "chr_" + Chr
                cmd = cmd + " " + chrtag + "/" + "chr" + Chr + ".PC.coordinates.txt"
            if results.mcomp is not None:
                cmd = cmd + " 2 "
            else:
                cmd = cmd + " 1 "
            toprint = "DIFF. COMPT. CALLING RUN ON ALL FULL CHROMOSOMES: "
            with open("chr_unaffected.txt", "w") as ucfile:
                for g in goodlist:
                    ucfile.write(g + "\n")
                    toprint = toprint + g + " "
            print(toprint)
            cmd = cmd + " chr_unaffected.txt " + str(results.res) + " " + distfile + " samplefile.txt " + scriptdir
            print(cmd)
            os.system(cmd)
            
            badlist = relist[1]
            for Chr in badlist:
                checkInputs(Chr)
                if results.reps is None:
                    getRepParams(Chr)
                makeSampleFile()
                chrtag = "chr_" + Chr
                cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -chr " + Chr + " -combined 1 -res " + results.res
                if results.blacklist is not None:
                    cmd = cmd + " -blacklist " + results.blacklist
                print(cmd)
                os.system(cmd)
                cmd = "Rscript " + os.path.join(scriptdir, "diffcmp_pythonV.r") 
                cmd = cmd + " " + chrtag + "/" + "chr" + Chr + ".PC.coordinates.txt"
                if results.mcomp is not None:
                    cmd = cmd + " 2 "
                else:
                    cmd = cmd + " 1 " 
                cmd = cmd + " " + Chr + " " + str(results.res) + " " + distfile + " samplefile.txt " + scriptdir
                print(cmd)
                os.system(cmd)  
            

if __name__ == '__main__':
    main()
