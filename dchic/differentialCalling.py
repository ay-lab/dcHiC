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

parser.add_argument("-inputFile", action = 'store', dest = 'inputfile', help = "input text file (right now is input.txt)")

parser.add_argument("-chrFile", action = 'store', dest = 'chr', help = "Chr file")

parser.add_argument("-makePlots", action = 'store', dest = 'makePlots', help = "Set to true if you want plots to be made. DEBUG feature. Set true by default.")

parser.add_argument("-res", action = 'store', dest = 'res', help = "Resolution")

parser.add_argument("-multiComp", action = 'store', dest = 'mcomp', help = "Include this field (-multiComp 1) if you want a combined comparison. DEBUG feature. Set true by default.")

#parser.add_argument("-keepIntermediates", action = 'store', dest = "keepIntermediates", help = "DEBUG FEATURE: Activate to output replicate fitting info while processing")

results = parser.parse_args()

startdir = os.getcwd()
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))

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

chrlist = []
with open(results.chr, "r") as input:
    for line in input:
        a = line.strip()
        chrlist.append(a)

if isGrouping == False:
    if (len(names) != len(groups)):
        print("Erroneous input file. Length of names, groups don't match.")
else:
    if (len(names) != len(groups)) or (len(groups) != len(groupings)):
        print("Erroneous input file. Length of names, groups, and groupings don't match.")
        
def setDestinations(Chr):
    global destinations
    destinations = []
    currdir = os.getcwd()
    chrTag = "chr_" + Chr
    pcFiles = os.path.join(currdir, chrTag, "pcFiles")
    pc_decision_file = os.path.join(currdir, chrTag, "pc_decision.txt")
    with open(pc_decision_file, "r") as pc_file:
        pc_decision = int(pc_file.readline().strip())
    if pc_decision == 1:
        for a in range(len(names)):
            filename = "pc1_" + str(names[a]) + "_exp_" + str(a+1) + ".txt"
            destinations.append(os.path.join(pcFiles, filename))
    if pc_decision == 2:
        for a in range(len(names)):
            filename = "pc2_" + str(names[a]) + "_exp_" + str(a+1) + ".txt"
            destinations.append(os.path.join(pcFiles, filename))

def checkInputs():
    global isGrouping
    with open(results.inputfile, 'r') as input:
        for line in input:
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
                    groupings_excl.append(temp[2])
                    groupings_sizes.append(0)
                groupings_sizes[groupings_excl.index(temp[2])] += 1
            else:
                print("Error in input file formatting. Exiting.")
                sys.exit(1)

def makeSampleFile():
    global names
    global groups
    print("Sample File in Progress")
    print(names)
    print(groups)
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
        plt.show()
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

def getRepParams(): 
    if results.makePlots is not None:
        makeDirectory()
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
                #print(len(r1data))
                #print(r1data[:10])
                #print(len(r2data))
                #print(r2data[:10])
                vals = createModel(r1data, r2data, str(Chr), False)
            for inc in range(len(r1data)):
                dist = (pointLineDist(r1data[inc], r2data[inc], vals[0], vals[1]))
                distances.append(dist)
            chrMean = statistics.mean(distances)
            chrStd = stdev(distances)
            print("Chromosome " + Chr + " /// Mean: " + str(chrMean) + " Std: " + str(chrStd))
            eucDistanceFile.write(str(chrMean) + "\t" + str(chrStd) + "\t" + str(Chr) + "\n")
    

def main():
    checkInputs()
    getRepParams()
    makeSampleFile()
    for Chr in chrlist:
        chrtag = "chr_" + Chr
        cmd = "python " + os.path.join(scriptdir, "makeBedGraph.py") + " -chr " + Chr + " -combined 1"
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
    cmd = cmd + " " + results.chr + " " + str(results.res) + " chrdistances.txt samplefile.txt " + scriptdir
    print(cmd)
    os.system(cmd)  

if __name__ == '__main__':
    main()
                   
                        
            