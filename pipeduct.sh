#!/bin/bash

######################################
#
#	usage pipeduct.sh project.txt fastq_dir
#
# project.txt is the project sample sheet
# fastqdir is the name of the folder where the fastq files are located
#
# assumes scripts are all on PATH (OK in pipeduct docker container
#	
#
#########################################


project='basename $1 ".txt"' #sample sheet for the project
fastqdir=$2



#########################################
#
# running all alignments
#
#########################################
bamdir=$project"_bam"

if [!-d $bamdir]; then
    mkdir -p $bamdir;
fi;

while read sampleID	replicateID	libraryID	filename	treatment_flag; 
do pipeduct_align.sh $fastqdir/$filename".fastq.gz" $bamdir
;done < $project

# get all mapping stats - generate Plot
grep '' $bamdir/*.stat | XXX	> $project.al -  #generate a table with filename total/mapped Nreads
alignQC.R 


#########################################
#
# running all bam2cnt
#
#########################################
cntdir=$project"_cnt"

if [!-d $cntdir]; then
    mkdir -p $cntdir;
fi;


while read sampleID	replicateID	libraryID	filename	treatment_flag; 
do pipeduct_count.pl $bamdir/$filename $bamdir
;done < $project


#########################################
#
# Import in R for plot and stats
#
#########################################
anadir=$project"_analysis"

if [!-d $anadir]; then
    mkdir -p $anadir;
fi;

pipeduct_analyze.R $project".txt" #assumes project.txt sample sheet and project_cnt folder for counts. 







