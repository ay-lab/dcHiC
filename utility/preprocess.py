#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:38:23 2020

@author: jeffreywang
"""

import argparse
import os
import sys
import hicstraw

parser = argparse.ArgumentParser()

parser.add_argument("-input", action= 'store', dest = 'input', help= "[required] 'cool' for .cool files, and 'hic' for .hic files.")

parser.add_argument("-file", action = "store", dest = "file", help = "[required] Enter the .cool/.mcool/.hic file path")

parser.add_argument("-res", action = "store", dest = "res", help = "[required] Enter the resolution for the processing")

parser.add_argument("-prefix", action = "store", dest = "prefix", help = "[required] Enter prefix of results")

parser.add_argument("-genomeFile", action = "store", dest = "genomeSize", help = "[.cool only] Location of chromosome size file (Must be: hg38, hg19, mm10, mm9)")

parser.add_argument("-removeChr", action = 'store', dest = 'remove', help = "[.hic only] Remove chromosomes: {Chr,Chr,Chr} format. Commonly used for Y. Default MT removed.")
            
results = parser.parse_args()

scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
res = int(results.res)

# .cool processing 

if results.input == "cool":
        
    # Extract chromosome size info

    chrs = []
    sizes = []
    chrSizes = {}
    with open(results.genomeSize, "r") as chrinput: # May need to change this to "rb" option if getting error 
        for line in chrinput:
            chrTag = line.strip().split()[0]
            chrSize = int(line.strip().split()[1])
            chrSizes[chrTag] = chrSize
            chrs.append(chrTag)
            sizes.append(chrSize)

    matrixArr = []
    filename, file_extension = os.path.splitext(results.file)

    # Set output names

    outname = results.prefix + "_" + results.res + ".matrix"
    extension = os.path.splitext(results.file)[1]
    filename = ""
    tmpfile = "tmp"
    if extension == ".mcool":
        filename = results.file + "::/resolutions/" + results.res 
    elif extension == ".cool":
        filename = results.file
    else:
        print("Extension is neither .cool or .mcool. Please make sure it is one of these two.")

    # Run cooler dump

    print("Creating sparse matrix............")
    cmd = "cooler dump -o " + tmpfile + " " + filename # cooler CLI one-index features don't work
    print(cmd)
    os.system(cmd)

    # Make one-indexed 

    with open(tmpfile, "r") as tmpIn:
        with open(outname, "w") as fileOut:
            for line in tmpIn:
                l = line.strip().split()
                fileOut.write(str(int(l[0])+1) + "\t" + str(int(l[1])+1) + "\t" + l[2] + "\n")

    name = results.prefix + "_" + str(results.res) + "_abs.bed" #data_200000_abs.bed
    iterator = 1 # one-based
        
    # Make bed files

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
        
else:

    # Extract Data

    hic = hicstraw.HiCFile(results.file)
    genome = hic.getGenomeID()
    resolutions = hic.getResolutions()

    print(f" - These are the resolutions of your file: {resolutions}")
    print(f" - This is the genome of your file: {genome}.")

    chrs = []
    sizes = []
    chrSizes = {}

    # Removing Chromosomes 

    chrsToRemove = []
    if results.remove:
        chrsToRemove = results.remove.strip().split(",") 
    chrsToRemove.append("all") # default
    chrsToRemove.append("mt")

    print(" - Removing these chromosomes:")
    for chrom in hic.getChromosomes():
        if chrom.name.lower() in chrsToRemove:
            if chrom.name.lower() == "all":
                print(f"  - {chrom.name} (.hic file artifact removed by default)")
            else:
                print(f"  - {chrom.name}")
        else:
            chrTag = "chr" + chrom.name
            chrSize = chrom.length
            chrSizes[chrTag] = chrSize
            chrs.append(chrTag)
            sizes.append(chrSize)

    # Make bed files

    name = results.prefix + "_" + str(results.res) + "_abs.bed" #data_200000_abs.bed
    iterator = 1 # one-based
    res = int(results.res)

    positionHash = {} # list of hash tables (first coord : index), one per chromosome

    with open(name, "w") as bedfile:
        print(" - Creating bed file")
        for a in range(len(chrs)):
            chrDic = {}
            length = chrSizes.get(chrs[a])
            posIterator = int(res)
            while posIterator <= length:
                chrDic[(posIterator-res)] = iterator
                bedfile.write(chrs[a] + "\t" + (str(posIterator-res)) + "\t" + str(posIterator) + "\t" + str(iterator) + "\n")
                posIterator+=res
                iterator+=1
            if (posIterator-res) < length:
                chrDic[(posIterator-res)] = iterator
                bedfile.write(chrs[a] + "\t" + (str(posIterator-res)) + "\t" + str(length) + "\t" + str(iterator) + "\n")
                iterator +=1
            positionHash[chrs[a]] = chrDic
        print(" - Bed File Creation Done.")

    # Sparse Matrix Dump

    outname = results.prefix + "_" + results.res + ".matrix"

    chrList = []
    print(" - Processing These Chromosomes: ")
    for c in range(len(chrs)):
        chr = chrs[c]
        chrSpecific = []
        chrNum = chr.split("chr")[1]
        print(f"  - {chr}")
        result = hicstraw.straw("observed", 'NONE', results.file, chrNum, chrNum, 'BP', res)
        for elem in result:
            chrSpecific.append(elem)
        chrSpecific.sort(key = lambda x : (x.binX, x.binY))
        chrList.append(chrSpecific)

    with open(outname, "w") as outfile:
        for c in range(len(chrList)):
            chr = chrs[c]
            for elem in chrList[c]:
                outfile.write("{0}\t{1}\t{2}\n".format(positionHash[chr][elem.binX], positionHash[chr][elem.binY], int(elem.counts)))
