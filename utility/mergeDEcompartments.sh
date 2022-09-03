#!/bin/bash
while getopts f:d: flag
do
    case "${flag}" in
        f) file=${OPTARG};;
        d) dist=${OPTARG};;
    esac
done

awk -v OFS='\t' '$(NF-1) < 0.1 {print $1,$2,$3}' $file |bedtools sort -i - |bedtools merge -d $dist -c 1 -o count |awk -v OFS='\t' '$NF > 1 {print $1,$2,$3}'
