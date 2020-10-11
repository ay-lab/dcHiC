#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 11:44:50 2020

@author: jeffreywang
"""

import argparse
import glob
import os
import math

parser = argparse.ArgumentParser()

parser.add_argument("-eigfile", action = 'store', dest = 'eigfile', help = "eigenvector file")

parser.add_argument("-chr", action = 'store', dest = 'chr', help = "chromosome num")

parser.add_argument("-exp", action = 'store', dest = 'exp', help = "experiment label")

parser.add_argument("-PC2", action = 'store', dest = 'pc2', help = "DEBUG FEATURE: Use this for PC2; different naming")

parser.add_argument("-combined", action = 'store', dest = 'combined', help = "Used for combination to large bedGraph in visualization step")

parser.add_argument("-res", action = 'store', dest = 'res', help = "Resolution.")

parser.add_argument("-blacklist", action = 'store', dest = 'blacklist', help = "Blacklist, if used.")

results = parser.parse_args()

if results.combined is None:
    outputLabel = "chr_" + results.chr + "_exp_" + results.exp
    chrLabel = "chr" + results.chr
    outputname = "HMFA_" + outputLabel + ".bedGraph"
    print(outputname)
    if results.pc2 is not None:
        outputname = "PC2_" + outputLabel + ".bedGraph"
    
    with open(results.eigfile, "r") as file:
        with open(outputname, "w") as output:
            trash = file.readline()
            #output.write("track name=\"" + outputLabel + "\"  yLineMark=\"0.0\" " + "alwaysZero=on maxHeightPixels=100:75:11 visibility=full viewLimits=-1:1 " + "autoScale=on type=bedGraph\n")
            line = file.readline()
            linez = line.split()
            lastBin = str(linez[0])[linez[0].index(".")+1:len(linez[0])]
            prevEig = linez[1] # this was to shift everything by one bin
            if int(lastBin) != 0:
                output.write(chrLabel + "\t0\t" + lastBin + "\t0\n")
            line = file.readline()
            while line:
                linez = line.split()
                currBin = str(linez[0])[linez[0].index(".")+1:len(linez[0])]
                output.write(chrLabel + "\t" + lastBin + "\t" + currBin + "\t" + str(prevEig) + "\n")
                prevEig = linez[1]
                lastBin = currBin
                line = file.readline()

else:
    if results.blacklist is not None:
        expRes = int(results.res)
        blacklistedregions = [] # only for this chr
        thisChrTag = "chr" + results.chr
        with open(results.blacklist, "r") as bfile:
            for line in bfile:
                l = line.strip().split()
                chrTag = l[0]
                if chrTag != thisChrTag:
                    continue
                start = int(l[1])
                end = int(l[2])
                badstart = math.ceil(start/expRes) * expRes
                badend = math.floor(end/expRes) * expRes
                for a in range(badstart, badend, expRes):
                    badbin = chrTag + "." + str(a)
                    blacklistedregions.append(badbin)
    # go into pcFiles
    chrtag = "chr_" + results.chr
    pcFiles = os.path.join(chrtag, "pcFiles")
    #pc_decision_file = os.path.join(chrtag, "pc_decision.txt")
    #with open(pc_decision_file, "r") as pc_file:
    #    pc_decision = int(pc_file.readline().strip())
    os.chdir(pcFiles)
    prefixes = []
    vals = []
    bins = []
    #if pc_decision == 1:
    #    globterm = "pc1"
    #else:
    #    globterm = "pc2"
    pcfilelist = glob.glob("pc_*")
    pcfilelist.sort(key=os.path.getmtime) # may be changed in the future 
    firstfile = True
    for file in pcfilelist:
        namearr = file.split("_")
        newstr = namearr[1]
        if "-" in namearr[1]:
            newstr = namearr[1].replace("-", ".")
        prefixes.append(newstr) # get the prefix
        with open(file, "r") as bg:
            bg.readline()
            bgvals = []
            for line in bg:
                line = line.strip().split()
                if results.blacklist is not None:
                    if line[0] in blacklistedregions:
                        continue
                bgvals.append(line[1])
                if firstfile:
                    bins.append(line[0])
                else:
                    continue
        if firstfile:
            firstfile = False
        vals.append(bgvals)
    os.chdir("..")
    combname = "chr" + results.chr + ".PC.coordinates.txt"
    #print((bins))
    for binpos in range(len(bins)):
        bins[binpos] = bins[binpos].replace(".", "-")
   # print(len(vals))
   # print(len(vals[1]))
   # print(len(prefixes))
   # print(len(bins))
    with open(combname, "w") as finalfile:
        for p in prefixes:
            finalfile.write(p + "\t")
        finalfile.write("\n")
        for a in range(len(bins)):
            finalfile.write(str(bins[a]) + "\t")
            for b in range(len(vals)):
                finalfile.write(vals[b][a] + "\t")
            finalfile.write("\n")
    os.chdir("..")
    
