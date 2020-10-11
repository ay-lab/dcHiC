#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 16:08:33 2020

@author: jeffreywang
"""
import argparse
import os
import copy
import re
import sys
import bisect
import glob
import math

parser = argparse.ArgumentParser() # 

parser.add_argument('-nExp', action = 'store', dest='totalExp', help='Number of experiments given')

parser.add_argument('-prePath', action='append', dest='tag_list',default=[], help='Give path to HOMER pre-processed temp file -prePath /to/exp1 -prePath /to/exp2 etc.')

parser.add_argument('-res', action = 'store', dest = 'res', help = 'resolution (10000, 25000, 50000, etc.)')

parser.add_argument("-numGroups", action = 'store', dest = 'numGroups', help = "number of groups entered")

parser.add_argument("-group", action = 'append', dest = 'groups', help = "number of experiments in a group. The sum of these should equal the number of experiments given. For instance: 8 experiments given, 3 groups created, group #s: 2, 3, 3")
                    
parser.add_argument("-expNames", action = 'append', dest = 'expNames', help = "a name for every experiment")

parser.add_argument("-groupNames", action = 'append', dest = 'groupNames', help = "a name for every group.")

parser.add_argument("-groupings", action = 'append', dest = 'groupings', help = "groupings")

parser.add_argument("-groupingNums", action = 'append', dest = 'groupingNums', help = "number of groups per grouping")

parser.add_argument('-chrNum', action = 'store', dest = 'chrNum', help = 'chromosome number: 1, 2, 3...X, Y')

parser.add_argument('-genome', action = 'store', dest = 'genome', help = "Used to determine eigenvector sign. Options: hg38, hg19, mm10, mm9 [If not provided, analysis will skip this step]")

#parser.add_argument('-signAnalysis', action = 'store', dest = 'analysis', help = "Used to determine which data will be used for determining eigenvector sign. Options: TSS, GC content. [If not provided, analysis will skip this step.]")

parser.add_argument("-alignData", action = 'store', dest = 'goldenpath', help = "Absolute location of GoldenPath data")

parser.add_argument('-grouping', action = 'store', dest = 'grouping', help = 'DEBUG FEATURE: Do average of eigenvectors (grouping) for comparisons -- 1 for true and 2 for false')

parser.add_argument("-blacklist", action = "store", dest = "blacklist", help = "ENCODE blacklisted regions - sorted, > 1mb")

parser.add_argument("-ncp", action = 'store', dest = "ncp", help = "ncp to use / scan in HMFA.")

# =============================================================================
# SYSTEM CHECKS
# =============================================================================

results = parser.parse_args()
numExp = int(results.totalExp)
numChr = 1
expRes = int(results.res)

if len(results.groups) != int(results.numGroups):
    print("INPUT ERROR: Number of 'group's != -numGroups")
    sys.exit()
if len(results.expNames) != numExp:
    print("INPUT ERROR: -nExp != Number of -expNames")
    sys.exit()
if len(results.groupNames) != int(results.numGroups):
    print("INPUT ERROR: -numGroups != Number of -groupNames")
    sys.exit()
if (int(results.numGroups) <= 1):
    print("Too few groups. EXIT")
    sys.exit()
#if len(results.groupingNums) != len(results.groupings):
#    print("Grouping Nums != # of Group Names Inputted.")
#    sys.exit()

sumOfGroupNums = 0
for a in results.groups:
    sumOfGroupNums += int(a)

if sumOfGroupNums != numExp:
    print("Group entering incorrect. Sum of -group's != -nExp")
    sys.exit()
    
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))

# =============================================================================
# RUNNING MFA
# =============================================================================

class chr_bin(object):
    pos = 0;
    def __init__(self, pos):
        self.pos = pos

def make_bin(pos):
    bin = chr_bin(pos)
    return bin

def IntersecOfSets(listoflists): 
    listoflistsbin = []
    for a in listoflists:
        newlist = []
        for x in a:
            try:
                newthing = make_bin(int(x[6:]))
                newlist.append(newthing)
            except:
                print("DEBUG:Blank space skipped. No worries.") # there are blank spaces within data at times, which are handled here. As stated, it is not a major problem. 
        listoflistsbin.append(newlist)
    newlistoflists = []
    for a in listoflistsbin:
        newlistoflists.append(set(a))
    finalset = [item for item in listoflists[0] if item in listoflists[1]]   
    i=0
    while i < (len(listoflistsbin)-1) :
        finalset = [item for item in finalset if item in listoflists[i+1]]
        i=i+1
    return finalset

listoflists = [] # [ [ [value, value, value], [value, value, value] ] ] <- there is an extra list "wrapper" because originally intended for multiple chromosomes. 
for i in range(numChr) :
    itemp = i + 1
    lists = []
    for j in range(numExp) :
        list = []
        jtemp = j + 1
        itempos = itemp * jtemp - 1
        obj = open(results.tag_list[itempos], "r")
        list = re.split(r'\t+', obj.readline())
        list.pop(0)
        a = copy.deepcopy(list)
        lists.append(a)
        obj.close()
    b = copy.deepcopy(lists)
    listoflists.append(b)

for a in listoflists: # handling new line characters in data
    for b in a:
        if len(b) == 0:
            print("Error with empty file. No values. Exiting.")
            sys.exit(1)
        b[len(b)-1] = b[len(b)-1].replace('\n', '')

commonlists = [] # one per chromosome
for i in listoflists:
    commonchrlist = IntersecOfSets(i)
    commonchrlist[len(commonchrlist)-1] = commonchrlist[len(commonchrlist)-1].replace('\n', '')
    commonlists.append(commonchrlist)

print("Positions Not Shared Across All Data Sets, Chromosome " + results.chrNum + ":\n")
listdiff = [] # All the positions not shared across data sets (array positions)
DifferentPositions= [] 
somethingIsDifferent = False
for l in listoflists:
    for a in l:
        tempdifflist = []
        for x in a:
            if x not in commonlists[0]: 
                #print(x)
                somethingIsDifferent = True
                if x not in DifferentPositions:
                    DifferentPositions.append(x)
                tempdifflist.append(a.index(x))
        listdiff.append(tempdifflist)
nameslist = []

DifferentPositions.sort()

if somethingIsDifferent == False:
    print("None")

else:
    for lineNum in DifferentPositions:
        print(lineNum)

for i in range(numChr):
    itemp = i + 1
    for j in range(numExp):
       jtemp = j + 1
       itempos = itemp * jtemp - 1
       newname = "BalancedChrMatrix_" + "exp_" + results.expNames[j] + ".txt" # this is done in the order originally presented via command line input
       nameslist.append(newname)
       with open(results.tag_list[itempos], "r") as input:
           with open(newname, "w") as output: 
               for element in commonlists[i]:
                   output.write("\t" + element)
               output.write("\n") 
               x = True
               if not listdiff[itempos]:
                   for line in input:
                       if x:
                           x= False
                       else:
                           output.write(line)
               else:
                   for line in input:
                       line = line.rstrip() # change made to get rid of \n at end # takes advantage of the fact that the first line has tab 
                       temp = line.split("\t")
                       if temp[0] in commonlists[i]: 
                           for x in range(len(temp)):
                               y=x-1
                               if y not in listdiff[itempos]:
                                   if x == 0:
                                       output.write(temp[x])
                                   else:
                                       output.write("\t" + temp[x])
                           output.write("\n") # changed to prevent eol error

removedpositions = []
if results.blacklist is not None:
    blacklistedregions = [] # only for this chr
    thisChrTag = "chr" + results.chrNum
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
                badbin = chrTag + "-" + str(a)
                blacklistedregions.append(badbin)
    #print(blacklistedregions)
    for i in range(numChr):
        itemp = i + 1
        for j in range(numExp):
           badregions = []
           badcoords = []
           jtemp = j + 1
           itempos = itemp * jtemp - 1
           lines = []
           name = "BalancedChrMatrix_" + "exp_" + results.expNames[j] + ".txt" 
           outline = ""
           with open(name, "r") as infile:
               line1 = infile.readline()
               l1data = line1.strip().split()
               for a in range(len(l1data)):
                   if l1data[a] in blacklistedregions:
                       badregions.append(a) # index
                       badcoords.append(l1data[a])
                   else:
                       outline = outline + "\t" + l1data[a]
               for line in infile:
                   lines.append(line.strip().split()[1:])
           removedpositions.append(badcoords)
           #print(badregions) 
           with open(name, "w") as newfile:
               newfile.write(outline + "\n")
               for a in range(len(lines)):
                   if a in badregions:
                       continue
                   else:
                       outline = l1data[a]
                       for b in range(len(lines[a])):
                           if b in badregions:
                               continue
                           else:
                               outline = outline + "\t" + lines[a][b]
                       newfile.write(outline + "\n") 
#print(removedpositions)
for a in range(len(commonlists)):
    commonlists[a] = [x for x in commonlists[a] if x not in removedpositions[a]]

commonbins = []
for bin_pos in commonlists[0]:
    commonbins.append(int(bin_pos[bin_pos.index("-")+1:]))
    
print(len(commonlists[a]))
cmd = "Rscript --vanilla " + (os.path.join(scriptdir, "Hier.R")) + " "

# Rscript --vanilla Hier.R BalancedChrMatrix_exp_1.txt BalancedChrMatrix_exp_2.txt BalancedChrMatrix_exp_3.txt BalancedChrMatrix_exp_4.txt 3 1 mcf7 mcf10a t47d hmec cancer normal 2 4 337

for tempvar1 in range(numExp):
    cmd = cmd + " " + nameslist.pop(0)
for tempvar2 in results.groups:
    cmd = cmd + " " + tempvar2
for tempvar3 in results.expNames:
    cmd = cmd + " " + tempvar3
for tempvar4 in results.groupNames:
    cmd = cmd + " " + tempvar4
cmd = cmd + " " + str(results.numGroups)
cmd = cmd + " " + str(numExp)
cmd = cmd + " " + str(len(commonlists[0]))
chrTag = "chr" + str(results.chrNum)

#if "p" in chrTag:
#    chrTag = chrTag[:chrTag.index("p")]
#if "q" in chrTag:
#    chrTag = chrTag[:chrTag.index("q")]

if results.genome is not None: 
    cmd = cmd + " " + str(results.res) + " " + chrTag + " " + str(results.genome) 
else:
    cmd = cmd + " None None None"
if results.goldenpath is None:
    cmd = cmd + " None"
else:
    cmd = cmd + " " + str(results.goldenpath)
if results.ncp is not None:
    cmd = cmd + " " + results.ncp
else:
    cmd = cmd + " " + 2
cmd = cmd + " " + scriptdir

print("\n" + cmd)
print("\nR OUTPUT:\n")
os.system(cmd)
print("\nR function done!\n")

oneChrList = []
for tv1 in range(numChr):
    chr_list = []
    for tv2 in range(numExp):
        name = "hmfa_chrRAW_" + results.expNames[tv2] + "_exp_" + str(tv2+1) + ".txt"
        dmfa_file = open(name, "r")
        dmfa_file.readline()
        exp_list = []
        for tv3 in range(len(commonlists[tv1])):
            el_list = dmfa_file.readline().split()
            exp_list.append(el_list[1])
        chr_list.append(exp_list)
    oneChrList = chr_list

chrExpBins = []

for topone in oneChrList:
    temp_list = []
    for bottomone in topone:
        temp_list.append(float(bottomone))
    chrExpBins.append(temp_list)
    
# =============================================================================
# Earlier, the balanced chromosome matrix is created by only taking common values. Those values not shared by all now need to be restored. This restores them by taking the average of the nearest existing values to any missing one.  
# =============================================================================
      
commonMax = max(commonbins)
commonMin = min(commonbins)
#completeCoverage = []
#PC1  
print(len(chrExpBins))
for ls in range(numExp):
    print("Removed Positions")
    print(removedpositions[ls])
    newName = "pc_" + results.expNames[ls] + "_exp_" + str(ls+1) + ".txt"
    with open(newName, "w") as hmfaFile:
        hmfaFile.write("\t\t\tDim.1\n")
        for i in range(commonMin, commonMax, expRes):
            firstColVal = "chr" + str(results.chrNum) + "." + str(i)
            otherFirstCol = "chr" + str(results.chrNum) + "-" + str(i)
            if otherFirstCol in removedpositions[ls]:
                #print("Removed " + otherFirstCol)
                hmfaFile.write(firstColVal + "\t0\n")
                continue
            if i in commonbins:
                pos = bisect.bisect_left(commonbins, i)
                hmfaFile.write(firstColVal + "\t" + str(chrExpBins[ls][pos]) + "\n")
            else:
                leftVal = bisect.bisect_left(commonbins, i) - 1
                rightVal = bisect.bisect_left(commonbins, i) 
                avgVal = (chrExpBins[ls][rightVal] + chrExpBins[ls][leftVal]) / 2
                hmfaFile.write(firstColVal + "\t" + str(avgVal) + "\n")  
for vt in range(numExp):
    name = "hmfa_chrRAW_" + results.expNames[vt] + "_exp_" + str(vt+1) + ".txt"
    cmd = "rm " + name
    os.system(cmd)
        
if int(results.grouping) == 2: # Debug Feature
    doNothing = True
    
else:
    oneChrList = []
    numBins = 0
    for tempvar in range(numExp):
        newName = "pc_" + results.expNames[tempvar] + "_exp_" + str(tempvar+1) + ".txt"
        a = True
        exp_list = []
        with open(newName, "r") as input:
            for line in input:
                if a:
                    a = False
                    continue
                else:
                    numBins += 1
                    el_list = line.split()
                    exp_list.append(el_list[1])
        oneChrList.append(exp_list)
        
    allBins = []
    newName = "pc_" + results.expNames[tempvar] + "_exp_" + str(tempvar+1) + ".txt"
    a = True
    with open(newName, "r") as input:
        for line in input:
            if a:
                a = False
                continue
            else:
                line1 = line.split()[0]
                tempchr = line1.split(".")[0]
                temploc = line1.split(".")[1]
                newtag = tempchr + "-" + temploc
                allBins.append(newtag)
    
    big_group_arr = []
    iterator = 0
    for num in results.groups: # this works by combining the eigenvectors
        size = int(num)
        groupArr = []
        b = True # first file in group
        for a in range(1, size+1):
            phrase = "pc_*"
            for file in glob.glob(phrase):
                filenum = int(file.split("_")[3].split(".")[0])
                if (a + iterator) == filenum:
                    tempvar = True
                    with open(file, "r") as input:
                        if b:
                            for line in input:
                                if tempvar:
                                    tempvar = False
                                    continue
                                else:
                                    groupArr.append(float(line.strip().split("\t")[1]))
                        else:
                            it = 0
                            for line in input:
                                if tempvar:
                                    tempvar = False
                                    continue
                                else:
                                    valToAdd = (float(line.strip().split("\t")[1]))
                                    groupArr[it] += valToAdd
                                    it += 1
                    if b:
                        b = False
        iterator += size
        for iterate in range(len(groupArr)):
            groupArr[iterate] /= float(size)
        big_group_arr.append(groupArr)
    
    for bigiterator in range(len(results.groups)): 
        newName = "hmfa_" + results.groupNames[bigiterator] + "_exp_" + str(bigiterator + 1) + ".txt"
        with open(newName, "w") as output:
            output.write("\t\t\tDim.1\n")
            print(len(big_group_arr))
            print(len(big_group_arr[0]))
            for i in range(commonMin, commonMax, expRes):
                firstColVal = "chr" + str(results.chrNum) + "." + str(i)
                output.write(firstColVal + "\t" + str(big_group_arr[bigiterator][int((i-commonMin)/expRes)]) + "\n")
    
    if results.groupings is not None:
        groupNum = 0                    
        for iterator in range(len(results.groupings)):
            newName = results.groupings[iterator] + "_grouping.txt"
            #numGroups = results.group
            with open(newName, "w") as output:
                output.write("\t\t\tDim.1\n")
                for i in range(commonMin, commonMax, expRes):
                    firstColVal = "chr" + str(results.chrNum) + "." + str(i)
                    totalVal = 0
                    for smalliterator in range(groupNum, groupNum + int(results.groupingNums[iterator])):
                        try:
                            totalVal += big_group_arr[smalliterator][int((i-commonMin)/expRes)]
                        except:
                            print(len(big_group_arr))
                            print(len(big_group_arr[0]))
                            print(smalliterator)
                            print(i)
                            print(commonMin)
                            print(expRes)
                    avgVal = float(totalVal) / int(results.groupingNums[iterator])
                    output.write(firstColVal + "\t" + str(avgVal) + "\n")
            groupNum += int(results.groupingNums[iterator])
    
    def isSameSign(a, b):
        a = float(a)
        b = float(b)
        state1 = a >= 0 and b >= 0
        state2 = a < 0 and b < 0
        return state1 or state2
    
    newChrList = []
    iterator = 0
    for i in range(len(results.groups)):
        groupList = []
        x = True
        for b in range(iterator, iterator + int(results.groups[i])):
            for c in range(len(allBins)): # changed to allBins, 2/25
                if x:
                    groupList.append(float(oneChrList[b][c]))
                else:
                    groupList[c] += (float(oneChrList[b][c])) 
            x = False
        for n in range(len(groupList)):
            groupList[n] /= float(results.groups[i])
        newChrList.append(groupList)
        iterator += int(results.groups[i])
    
    with open("lengthdoc.txt", "w") as out:
        out.write(str(len(allBins))+"\n")
        pos1 = str(allBins[0]).split("-")[1]
        pos2 = str(allBins[len(allBins)-1]).split("-")[1]
        out.write(str(pos1) + "\t" + str(pos2) + "\n")
        #out.write(str(pcNum)) # Used in visualize.py
    
    for a in range(len(results.groups)): ### NEED TO CHANGE BECAUSE USES newChrList / oneChrList
        for b in range(a+1, len(results.groups)):
            name = "compartmentSwitch_" + results.groupNames[a] + "_" + results.groupNames[b] + "_" + str(a+1) + "_" + str(b+1) + ".txt"
            with open(name, "w") as out:
                for x in range(len(newChrList[0])):
                    if not isSameSign(newChrList[a][x], newChrList[b][x]):
                        out.write(allBins[x] + "\n")
