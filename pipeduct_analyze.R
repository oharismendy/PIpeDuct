######################################
#
#	usage pipeduct_analyze.R sample_sheet.txt
#
# 
#########################################
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)

#args <- commandArgs(TRUE)
#project <- read.delim(args[1])
#projectID<-sub(".txt","",args[1])

args<-"mytestproject.txt"
projectID<-sub(".txt","",args)
project <- read.delim("~/Documents/UCSD/Admin/IT/PipeDuct/test_data/sample_sheet.txt")
project$cntfiles<-trimws(paste(getwd(),paste(sub(".txt","",args),"count",sep="_"),sub(".fastq.gz",".cnt",project$filename),sep="/"))

#import all files 
read.tables <- function(file.names) {
  require(plyr)
  ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn,header=FALSE,sep="\t")))
}

data <- read.tables(project$cntfiles)
colnames(data)<-c("cntfiles","locus","chr","coord","strand","cnt","seq")
data<-select(left_join(data,project),X.sampleID,replicateID,libraryID,treatment_flag,locus,cnt,seq)

######################################
#
# getting some QC and	count normalization and QC
#
# 
#########################################

data$Ncnt<-ave(data$cnt,data$libraryID,FUN=function(x) x/sum(x))


stat<-select(project,libraryID)
stat<-left_join(stat,aggregate(locus~libraryID,data=data,FUN=length))
stat<-left_join(stat,aggregate(cnt~libraryID,data=data,FUN=sum))
stat<-left_join(stat,dplyr::rename(aggregate(cnt~libraryID,data=data,function(x) length(x[x>1])/length(x)),frac_rec_loci=cnt))
stat<-left_join(stat,dplyr::rename(aggregate(Ncnt~libraryID,data=data,function(x) quantile(x,0.75)),top_quartile=Ncnt))
stat<-left_join(stat,dplyr::rename(aggregate(Ncnt~libraryID,data=data,function(x) quantile(x,1)),max=Ncnt))
stat<-left_join(stat,dplyr::rename(aggregate(Ncnt~libraryID,data=data,function(x) quantile(x,0.5)),median=Ncnt))

write.table(stat,paste(projectID,"cnt_stat","txt",sep="."),sep="\t",quote=FALSE,row.names = FALSE)



######################################
#
# overlaps between libraries
#
# 
#########################################

recur<-dcast(locus~libraryID,data=data,value.var="Ncnt",fun.aggregate=mean)
rownames(recur)<-recur$locus
recur<-select(recur,-locus)
recur<-filter(recur,rowSums(!is.na(recur))>1)




######################################
#
#	nucleotide context
#
# 
#########################################

