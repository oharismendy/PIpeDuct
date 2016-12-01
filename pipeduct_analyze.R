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

args <- commandArgs(TRUE)
project <- read.delim(args[1])
projectID<-sub(".txt","",basename(args[1]))
projectdir<-dirname(args[1])

project$cntfiles<-paste(projectdir,paste(projectID,"cnt",sep="_"),sub(".fastq.gz",".cnt",project$filename),sep="/")

#import all files 
read.tables <- function(file.names) {
  require(plyr)
  ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn,header=FALSE,sep="\t")))
}

data <- read.tables(project$cntfiles)
colnames(data)<-c("cntfiles","locus","chr","coord","strand","cnt","seq")
data<-left_join(data,project)

######################################
#
# getting some QC and	count normalization and QC
#
# 
#########################################

data$Ncnt<-ave(data$cnt,data$libraryID,FUN=function(x) x/sum(x))
data$tile<-ave(data$Ncnt,data$libraryID,FUN=function(x) ntile(x,10))

#adding sequence info
data<-filter(data,length(data$seq)>7)
data<-mutate(data,AG=ifelse(substr(data$seq,5,6)=="AG" | substr(data$seq,5,6)=="GG","[AG]G","not[AG]G"))
data<-mutate(data,dinuc=substr(data$seq,5,6))

stat<-select(project,libraryID)
stat<-left_join(stat,aggregate(locus~libraryID,data=data,FUN=length))
stat<-left_join(stat,aggregate(cnt~libraryID,data=data,FUN=sum))
stat<-left_join(stat,dplyr::rename(aggregate(cnt~libraryID,data=data,function(x) length(x[x>1])/length(x)),frac_rec_loci=cnt))
stat<-left_join(stat,dplyr::rename(aggregate(Ncnt~libraryID,data=data,function(x) quantile(x,0.75)),top_quartile=Ncnt))
stat<-left_join(stat,dplyr::rename(aggregate(Ncnt~libraryID,data=data,function(x) quantile(x,1)),max=Ncnt))
stat<-left_join(stat,dplyr::rename(aggregate(Ncnt~libraryID,data=data,function(x) quantile(x,0.5)),median=Ncnt))

write.table(stat,paste(projectID,"cnt_stat","txt",sep="."),sep="\t",quote=FALSE,row.names = FALSE)

pdf("cnt_dist.pdf")
ggplot(data,aes(log10(Ncnt),color=libraryID))+stat_ecdf()+xlab("Normalized Count (log10)")+ylab("Fraction of loci")
#dev.off()

######################################
#
# overlaps between libraries
#
# 
#########################################

#recur<-dcast(locus~libraryID,data=data,value.var="Ncnt",fun.aggregate=mean)
#rownames(recur)<-recur$locus
#recur<-select(recur,-locus)
#recur<-filter(recur,rowSums(!is.na(recur))>1)
#
#recur<-dcast(locus+AG~libraryID,data=data,value.var="Ncnt",fun.aggregate=mean)
#
#recur<-mutate(recur,sharingControl=ifelse(!is.na(cDDP_rep2) & !is.na(untreated),"shared","unique"))
#recur<-mutate(recur,sharingControl=ifelse(sharingControl=="shared" & cDDP_rep2>5*untreated,"cDDPhigh",sharingControl))
#recur<-mutate(recur,sharingCDDP=ifelse(!is.na(cDDP_rep2) & !is.na(cDDP_rep1) & is.na(untreated),"shared","unique"))
#recur<-mutate(recur,uniqueControl=ifelse(!is.na(untreated) & is.na(cDDP_rep1) & is.na(cDDP_rep2),"unique","shared"))

grid<-expand.grid(unique(data$libraryID),unique(data$libraryID))
colnames(grid)<-c("libraryID1","libraryID2")

overlap<-function(x,y) {
  return(nrow(inner_join(filter(data,libraryID==x),filter(data,libraryID==y),by="locus")))
}

overlapAG<-function(x,y) {
  return(nrow(inner_join(filter(data,AG=="[AG]G" & libraryID==x),filter(data,AG=="[AG]G" & libraryID==y),by="locus")))
}

grid$cnt<-mapply(overlap,as.character(grid[,1]),as.character(grid[,2]),SIMPLIFY = TRUE)
grid$cntAG<-mapply(overlapAG,as.character(grid[,1]),as.character(grid[,2]),SIMPLIFY = TRUE)
write.table(grid,paste(projectID,"overlap_stat","txt",sep="."),sep="\t",quote=FALSE,row.names = FALSE)


################################################
#
# nucleotide context by coverage quantile. 
# for each dinculeotide, calculate the Odds Ratio and p-value. at different N-cov cutoff. 
#
# 
###############################################

pdf("dinuc.pdf")
#dinuc_agg<-aggregate(locus~libraryID+tile+dinuc,data=data,FUN=length)
#write.table(dinuc_agg,paste(projectID,"dinuc_stat","txt",sep="."),sep="\t",quote=FALSE,row.names = FALSE)

dinuc_agg<-aggregate(locus~libraryID+treatment_flag+dinuc,data=filter(data,tile==10),FUN=length)
dinuc_agg$total<-ave(dinuc_agg$locus,dinuc_agg$libraryID,FUN=sum)
dinuc_agg$Frac<-dinuc_agg$locus/dinuc_agg$total
dinuc_agg<-filter(dinuc_agg,!grepl("N",dinuc_agg$dinuc))

ggplot(dinuc_agg,aes(dinuc,Frac,fill=dinuc))+
	geom_bar(stat="identity")+
	facet_wrap(~libraryID,ncol=1)


######################################
#
#	Export BED files for enrichment analysis
#
# 
#########################################



exportBED<-function(x){
	tmp<-filter(data,libraryID==x & tile==10 & AG=="[AG]G")	
	write.table(paste(tmp$chr,tmp$coord-1,tmp$coord,tmp$locus,"100",tmp$strand,sep="\t"),paste(x,"bed",sep="."),col.names=FALSE, quote=FALSE,row.names = FALSE,sep="\t")
	
}

mapply(exportBED,unique(data$libraryID),SIMPLIFY = TRUE)

dev.off()


######################################
#
#	ratio to control and precision
#
# 
#########################################



######################################
#
# investigate interstrands loci 
#
# 
#########################################
