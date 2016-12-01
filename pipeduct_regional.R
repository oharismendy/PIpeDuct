library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
project <- read.delim(args[1])
projectID<-sub(".txt","",basename(args[1]))
projectdir<-dirname(args[1])
reg <- args[2]
regID<-sub(".bed","",basename(args[2]))


#########################################
#
# import the intersection counts
#
#########################################

project$interfiles<-paste(projectdir,paste(projectID,"regional",sep="_"),paste(project$libraryID,regID,"txt",sep="."),sep="/")

read.tables <- function(file.names) {
  require(plyr)
  ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn,header=FALSE,sep="\t")))
}

data <- read.tables(project$interfiles)
colnames(data)<-c("interfiles","region","cnt")
data<-left_join(data,project)

#########################################
#
# sum region counts of all untreated libraries 
#
#########################################

dataUntreated<-aggregate(cnt~region, data=filter(data,treatment_flag==FALSE), FUN=sum)
dataUntreated$libraryID<-"untreated"

#########################################
#
# test each region x treated library 
#
#########################################
testregion<-function(x,y){
  
	t<-rbind(select(filter(data,libraryID==x),libraryID,region,cnt),dataUntreated)
	t<-mutate(t,category=ifelse(region==y,"in","out"))
	t_agg<-aggregate(cnt~category+libraryID,data=t,FUN=sum)
	
	pval<-chisq.test(matrix(c(filter(t_agg,libraryID==x & category=="in")$cnt,filter(t_agg,libraryID==x & category=="out")$cnt,filter(t_agg,libraryID=="untreated" & category=="in")$cnt,filter(t_agg,libraryID=="untreated" & category=="out")$cnt),2,2))$p.value
	ratio<-(filter(t_agg,libraryID==x & category=="in")$cnt/filter(t_agg,libraryID==x & category=="out")$cnt)/(filter(t_agg,libraryID=="untreated" & category=="in")$cnt/filter(t_agg,libraryID=="untreated" & category=="out")$cnt)	
	return(c(pval,ratio))
}
	
volcano<-expand.grid(unique(filter(data,treatment_flag=="TRUE")$libraryID),unique(data$region))
tmp<-mapply(testregion,as.character(grid[,1]),as.character(grid[,2]),SIMPLIFY = TRUE)
volcano<-cbind(volcano,t(tmp))
colnames(volcano)<-c("libraryID","region","pval","ratio")

write.table(volcano  "volcano.txt",quote=FALSE,row.names = FALSE,sep="\t")
	
#########################################
#
# volcano plot
#
#########################################

pdf("volcano.pdf")

ggplot(volcano,aes(log2(ratio),-log10(pval),color=region,shape=libraryID))+
  geom_point(size=5)+
  xlim(c(-1,1))+
  xlab("Treated / Untreated (log2)")+
  ylab("-log10(P.value)")+
  theme(text=element_text(size=26))
       
dev.off()

