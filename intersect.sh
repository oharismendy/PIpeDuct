#!/bin/bash

# this command uses Bedtools intersectbed to determine overlaps with regions in a bed file

$adduct_pos=$1 # a ".txt: file listing all addcut postions to be examined. in a chr,position format
$regions=$2 #BED formatted list of regions to examine chr,start,stop,region_name

name=`basename $adduct_pos ".txt"`

cat $file | awk '{start=$2-1; print $1,start,$2}' | sed 's/ /\t/g' | grep '^chr' | intersectBed -a $regions -b - -wo | cut -f 4 | sort | uniq -c | awk '{print $2,$1}' > $name.intersect.txt
