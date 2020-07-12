#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:38:23 2020

@author: jeffreywang
"""

import argparse
import os
import sys

parser = argparse.ArgumentParser()

parser.add_argument("-input", action = "store", dest = "file", help = "Enter path to pre-processed Pearsons matrix if -mode = 1. If -mode = 2, enter the .cool/.mcool file")

parser.add_argument("-mode", action = "store", dest = "mode", help = "Enter 1 (pre-processed Pearsons matrix from Juicebox) or 2 (.cool bed file)")

parser.add_argument("-genomeFile", action = "store", dest = "genomeSize", help = "Location of chromosome size file (Must be: hg38, hg19, mm10, mm9)")

parser.add_argument("-res", action = "store", dest = "res", help = "Enter the resolution for the processing")

parser.add_argument("-prefix", action = "store", dest = "prefix", help = "FOR OPTION 1: Enter chromosome of file (i,e, 'X' or '3'). FOR OPTION 2: Enter prefix of file you wish to run")
                    
results = parser.parse_args()

scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
res = int(results.res)

chrs = []
sizes = []
chrSizes = {}
with open(results.genomeSize, "r") as chrinput:
    for line in chrinput:
        chrTag = line.strip().split()[0]
        chrSize = int(line.strip().split()[1])
        chrSizes[chrTag] = chrSize
        chrs.append(chrTag)
        sizes.append(chrSize)

matrixArr = []

if results.mode == "1": #.hic file data dump
    
    chromosome = results.prefix
    if chromosome == "x" or chromosome == "y":
        chromosome = chromosome.upper()
    chrTag = "chr" + chromosome
    
    with open(results.file, "r") as locfile:
        for line in locfile:
            matrixArr.append(line.split())
            
    lineEls = []
    badPositions = [] # 0-indexed
    
    for line in matrixArr:        
        goodLine = False
        perfectLine = True
        for a in line:
            if a != "NaN":
                goodLine = True
            if a == "NaN":
                perfectLine = False
        if goodLine:
            for a in range(len(line)):
                if line[a] == "NaN":
                    badPositions.append(a)
            break
    
    newName = "chr" + str(results.chr) + ".matrix"
    iterator = 0
    outline = ""
    while iterator < len(matrixArr):
        if iterator in badPositions:
            iterator+=1
            continue
        outline += "\t" + chrTag + "-" + str(res*iterator)
        iterator+=1
        
    with open(newName, "w") as outfile:
        outfile.write(outline + "\n")
        for a in range(len(matrixArr)):
            if a in badPositions:
                outline = chrTag + "-" + str(res*a)
                for b in range(len(line)):
                    outline += "\t0"
            else:
                outline = chrTag + "-" + str(res*a)
                line = matrixArr[a]
                for b in range(len(line)):
                    if b in badPositions:
                        outline += "\t0"
                    else: 
                        outline += "\t" + str(line[b])
                outfile.write(outline + "\n")     

elif results.mode == "2": #.cool data dump option
    name = results.prefix + "_" + str(results.res) + "_abs.bed" #data_200000_abs.bed
    iterator = 0
        
    with open(name, "w") as bedfile:
        print("Creating bed file.............")
        for a in range(len(chrs)):
            length = chrSizes.get(chrs[a])
            posIterator = res
            while posIterator <= length:
                bedfile.write(chrs[a] + "\t" + (str(posIterator-res)) + "\t" + str(posIterator) + "\t" + str(iterator) + "\n")
                posIterator+=res
                iterator+=1
            if (posIterator-res) < length:
                bedfile.write(chrs[a] + "\t" + (str(posIterator-res)) + "\t" + str(length) + "\t" + str(iterator) + "\n")
                iterator +=1
    
    print("Creating sparse matrix............")
    outname = results.prefix + "_" + results.res + ".matrix"
    extension = os.path.splitext(results.file)[1]
    filename = ""
    if extension == ".mcool":
        filename = results.file + "::/resolutions/" + results.res 
    elif extension == ".cool":
        filename = results.file
    else:
        print("Extension is neither .cool or .mcool. Please make sure it is one of these two.")
    
    cmd = "cooler dump -o " + outname + " " + filename 
    os.system(cmd)
    
    
    
    
    
#    with open(results.file, "r") as matrixfile:
#        print("Editing sparse matrix............")
#        matrixfile.readline()
#        newname = "data_" + str(results.res) + ".matrix"
#        with open(newname, "w") as newmatrixfile:
#            for line in matrixfile:
#                linearr = line.split()
#                newmatrixfile.write(str(int(linearr[0])+1) + "\t" + str(int(linearr[1])+1) + "\t" + linearr[2] + "\n")
