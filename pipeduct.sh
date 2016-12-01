#!/bin/bash

######################################
#
#	usage pipeduct.sh project.txt fastq_dir
#
# project.txt is the project sample sheet
# fastqdir is the name of the folder where the fastq files are located
#
# assumes scripts are all on PATH (OK in pipeduct docker container)
#	
#
#########################################


project=`basename $1 ".txt"` #sample sheet for the project
projectdir="/mnt/oncogxA/Projects/TDTSEQ/NGS/PipeDuct_dev/"
#fastqdir=$2
scriptsdir="/mnt/oncogxA/Projects/TDTSEQ/NGS/PipeDuct_dev/scripts"


#########################################
#
# running all alignments
#
#########################################
bamdir=$projectdir/$project"_bam"

#if [[ ! -d $bamdir ]]; then
#    mkdir -p $bamdir
#fi
#
#
#while read sampleID	replicateID	libraryID	filename	treatment_flag; do 
#	$scriptsdir/pipeduct_align.sh $fastqdir/$libraryID".fastq.gz" $bamdir 2> $bamdir/$libraryID".align.log"
#done < $projectdir/$project".txt"

# get all mapping stats - generate Plot
#grep '' $bamdir/*.stat | XXX	> $project.al -  #generate a table with filename total/mapped Nreads
#alignQC.R 


#########################################
#
# running all bam2cnt
#
#########################################
cntdir=$projectdir/$project"_cnt"

if [[ ! -d $cntdir ]]; then
    mkdir -p $cntdir
fi

#assumues sudo pip install weblogo

#cd $cntdir
#while read sampleID	replicateID	libraryID	filename	treatment_flag; do 
#	if [[ ! $filename =~ filename ]] ;then 
#	$scriptsdir/pipeduct_count.pl $bamdir/$libraryID".bam" 2> $libraryID".count.log" & 
#	fi
#	done < $projectdir/$project".txt"

#wait

#while read sampleID	replicateID	libraryID	filename	treatment_flag; do 
#	cat $libraryID".cnt" | cut -f 1,6 | sed 's/chr/>chr/g' | sed 's/\t/\n/g' | weblogo -F pdf -i -4 -t $libraryID -a 'ACGT' -A 'dna' -c classic -S 0.1 > $libraryID".pdf"
#	cat $libraryID".cnt" | cut -f 1,6 | sed 's/chr/>chr/g' | sed 's/\t/\n/g' | weblogo -F pdf -i -4 -t $libraryID -a 'ACGT' -A 'dna' -c classic > $libraryID".full.pdf"
#done < $projectdir/$project".txt"


#cd $projectdir
#########################################
#
# getting Weblogo
#
#########################################
#http://weblogo.threeplusone.com/manual.html


#########################################
#
# Import in R for plot and stats
#
#########################################
anadir=$projectdir/$project"_analysis"

#if [[ ! -d $anadir ]]; then
#    mkdir -p $anadir
#fi
#
#cd $anadir
#Rscript $scriptsdir/pipeduct_analyze.R $projectdir/$project".txt" 2>$project".analyze.log" #assumes project.txt sample sheet and project_cnt folder for counts. 
#cd ..


#########################################
#
# Perfrom regional analysis 
#
#########################################
regdir=$projectdir/$project"_regional"

if [[ ! -d $regdir ]]; then
    mkdir -p $regdir
fi


regbed=hub_1623_HMM15IMR90_Cell_Line.bed
regname=`basename $regbed ".bed"`

while read sampleID	replicateID	libraryID	filename	treatment_flag; do 
bedtools sort -i $anadir/$libraryID.bed | intersectBed -wao -a - -b $projectdir/$regbed | cut -f 10 | sort | uniq -c | awk '{print $2"\t"$1}' > $regdir/$libraryID.$regname.txt
done < $projectdir/$project".txt"


Rscript $scriptsdir/pipeduct_regional.R $projectdir/$project".txt" $regbed 2>$project".regional.log"

