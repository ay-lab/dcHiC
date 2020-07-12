#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 11:43:57 2020

@author: jeffreywang
"""

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("-nExp", action = 'store', dest = 'numExp', help = "Number of experiments")

parser.add_argument("-nGroups", action = 'store', dest = 'numGroups', help = "Number of groups (if M vs M)")

parser.add_argument("-groups", action = "append", dest = "groupNames", help = "Names of groups")

parser.add_argument("-exp", action = 'append', dest = 'expNames', help = "append a name") # same order as exp1, exp2, etc.  -- no longer used

parser.add_argument("-chr", action = 'append', dest = 'chrs', help = "add chromosomes")

results = parser.parse_args()

nExp = int(results.numExp)

# =============================================================================
# Obtaining the coordinates from each chromosome's coordinate.txt
# =============================================================================

xcoords = [] # [[chr1: x1, x2, x3], [chr2: x1, x2, x3], etc.]
ycoords = []

for chrm in results.chrs:
    target_dir = "chr_" + chrm
    rel_path = os.path.join(target_dir, "coordinates.txt")
    xtemp = []
    ytemp = []
    with open(rel_path, "r") as input:
        a = True
        for line in input:
            arr_vals = line.split()
            if arr_vals[0] == "Dim.1" and not a:
                break
            if arr_vals[0] == "Dim.1":
                a = False
                continue
            else:
                xtemp.append(float(arr_vals[1]))
                ytemp.append(float(arr_vals[2]))
    xcoords.append(xtemp)
    ycoords.append(ytemp)

aggregate_xcoords = []
aggregate_ycoords = []

for a in range(int(results.numExp)):
    temp_xsum = 0
    temp_ysum = 0
    for b in range(len(xcoords)):
        temp_xsum += xcoords[b][a]
        temp_ysum += ycoords[b][a]
    aggregate_xcoords.append(temp_xsum)
    aggregate_ycoords.append(temp_ysum)

for a in range(int(results.numExp)):
    aggregate_xcoords[a] /= len(results.chrs)
    aggregate_ycoords[a] /= len(results.chrs)

# =============================================================================
# Creating figure (all experiments)
# =============================================================================

plt.figure(0)
plt.axis([0, 1, 0, 1])
plt.title("PC1 vs PC2, aggregate")
plt.xlabel("PC1")
plt.ylabel("PC2")

colors = []
for a in range(int(results.numExp)):
    rgb = np.random.rand(3,)
    colors.append(rgb)

for a in range(int(results.numExp)):
    plt.plot(aggregate_xcoords[a], aggregate_ycoords[a], 'ro', marker = "*", markersize = 15, color = colors[a])
    for b in range(len(xcoords)):
        plt.plot(xcoords[b][a], ycoords[b][a], marker = ".", markersize = 2, color = colors[a])

legend_els = []

for a in range(int(results.numExp)):
    legend_els.append(Line2D([0], [0], color=colors[a], lw=2, label=results.expNames[a]))

plt.legend(handles=legend_els, loc='upper right')
    
plt.savefig('experiment_coordinates.png', bbox_inches='tight', dpi = 600)
plt.close()

# =============================================================================
# Creating figure (group-wise)
# =============================================================================

group_x = []
group_y = []
group_aggx = []
group_aggy = []
for chrm in results.chrs:
    target_dir = "chr_" + chrm
    rel_path = os.path.join(target_dir, "coordinates.txt")
    xtemp = []
    ytemp = []
    with open(rel_path, "r") as input:
        a = True
        b = False
        for line in input:
            arr_vals = line.split()
            if arr_vals[0] == "Dim.1" and not a:
                b = True
            if arr_vals[0] == "Dim.1":
                a = False
                continue
            if b and arr_vals[0][:3] != "Dim":
                xtemp.append(float(arr_vals[1]))
                ytemp.append(float(arr_vals[2]))
        group_x.append(xtemp)
        group_y.append(ytemp)
    
#print(group_x)
#print(group_y)

for a in range(int(results.numGroups)):
    temp_xsum = 0
    temp_ysum = 0
    for b in range(len(group_x)):
        temp_xsum += group_x[b][a]
        temp_ysum += group_y[b][a]
    group_aggx.append(temp_xsum)
    group_aggy.append(temp_ysum)

for a in range(int(results.numGroups)):
    group_aggx[a] /= len(results.chrs)
    group_aggy[a] /= len(results.chrs)

for a in range(int(results.numGroups)):
    plt.figure(a+1)
    plt.title("PC1 vs PC2 for " + results.groupNames[a])
    plt.axis([0, 1, 0, 1])
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.plot(group_aggx[a], group_aggy[a], 'ro', marker = "*", markersize = 15, color = colors[a])
    for b in range(len(group_x)):
        plt.plot(group_x[b][a], group_y[b][a], marker = ".", markersize = 2, color = colors[a])
        if len(xcoords) > 1:
            plt.text(group_x[b][a] * (1+0.0035), group_y[b][a] * (1+0.0035), results.chrs[b], fontsize = 8)
    name = results.groupNames[a] + "_fig.png"
    plt.savefig(name, bbox_inches='tight', dpi = 600)
    plt.close()
    
    
print("Coordinate grids created.")
