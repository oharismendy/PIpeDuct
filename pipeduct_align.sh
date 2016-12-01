#!/bin/bash

######################################
#
#	usage align.sh file.fastq.gz outdir 
#
# parameters include 
#	reference genome fasta file
#
# dependencies include:
# 		bwa
#		samtools
#
#########################################


ref=/scratch/references/hg19/bwa_v7/hg19_lite.fa
name=`basename $1 ".fastq.gz"`
workdir=`dirname $1`
outdir=$2


bwa aln -t 3 $ref $workdir/$name.fastq.gz > $workdir/$name.sai
wait
bwa samse $ref $workdir/$name.sai $workdir/$name.fastq.gz | samtools view -buSh - > $workdir/$name.bam 
wait
samtools sort -m 1G $workdir/$name.bam $workdir/$name.sorted  
wait
mv $workdir/$name.sorted.bam $workdir/$name.bam
samtools flagstat $workdir/$name.bam > $workdir/$name.stat &
samtools index $workdir/$name.bam 

wait

rm $workdir/$name.sai
mv $workdir/$name.bam* $outdir
mv $workdir/$name.stat $outdir

 
