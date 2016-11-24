#!/usr/bin/perl -w

######################################
#
#	this script take a BAM file as input and outputs the locus_ID,chr,coord,strand,count,seq_context 
#
# parameters include 
#	MQ mapping quality
#	readlength
#	fastalen	number of nuceltodie context to display
#	reference genome fasta file
#
# dependencies include:
#	samtools
#	bedtools
#
#########################################


use strict; 
use File::Basename;


# set parameters
my $MQ=10 ;		#read miminal mapping quality
my $readlength=50;
my $fastalen=4; #number of nucleotides to grab on eitherside of the last read nculeotide. The total length will be therefore 2*fastlen+1
my $ref="/mnt/idash/Genomics/data_resources/references-and-indexes/hg19/hg19.fa"
my $name=basename($ARGV[0],".bam");

# declare variables
my ($start,$end,$strand);
my ($chr,$coord,$flag,$cnt);
my $bed;
my $fasta;

#########################################
#
# Get help
#
#########################################

if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
print "USAGE: bam2cnt.pl [file.bam] \n";
exit;
}

#########################################
#
# count reads per locus
#
#########################################

open (my $rawcount, "samtools view $ARGV[0] | awk '\$5>$MQ' | cut -f 2,3,4 | sort | uniq -c | awk '{print \$3,\$4,\$2,\$1}' | sed 's/ /\t/g' |");
open OUT, ">$name.cnt";


#########################################
#
# extract sequence context
#
#########################################
while (<$rawcount>){
chomp;
($chr,$coord,$flag,$cnt)=split;

if ($flag==0){
$start=$coord-$fastalen-1;
$end=$coord+$fastalen;
$strand="+";
}
else {
$start=$coord+$readlength-2-$fastalen;
$end=$coord+$readlength-1+$fastalen;
$strand="-";
}


$bed=$chr."\t".$start."\t".$end."\t"."strand\t1\t".$strand;
$fasta=`echo $bed | sed 's/ /\t/g' | bedtools getfasta -name -s -fi $ref -bed - -fo stdout | tail -n 1`;
$fasta=uc($fasta);

#########################################
#
# print out results
#
#########################################

chomp($fasta);
print OUT "$chr:$coord:$strand\t$chr\t$coord\t$strand\t$cnt\t$fasta\n";
}
close OUT;

