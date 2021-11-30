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

parser.add_argument("-input", action = "store", dest = "file", help = "Enter the .cool/.mcool file")

parser.add_argument("-genomeFile", action = "store", dest = "genomeSize", help = "Location of chromosome size file (Must be: hg38, hg19, mm10, mm9)")

parser.add_argument("-res", action = "store", dest = "res", help = "Enter the resolution for the processing")

parser.add_argument("-prefix", action = "store", dest = "prefix", help = "Enter prefix of results")
                    
results = parser.parse_args()

scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
res = int(results.res)

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
    
