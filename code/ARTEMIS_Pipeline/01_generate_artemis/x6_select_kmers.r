args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

library(data.table)
library(tidyverse)
library(devtools)
load_all("/useful.stuff.aa")           # Load

f<-list.files("./denovo_kmer_type",pattern=".txt",full.names=T)
sv<-gsub("./denovo_kmer_type","./unique_kmer_type",f[i])

if (!file.exists(sv)) {
	print("Work on it!")
	d<-fread(f[i])
	print(nrow(d))

	d<- d %>% rename("kmer"="V1")
	d<- d %>% rename("count"="V2")
	d$other<-0
	print("start!")
	f2<-f[-i]
	for (j in 1:length(f2)) {
  	dat<-squish(f2[j])
  	d<-left_join(d,dat,by=c("kmer"="V1"))
  	d[is.na(d)] <- 0  
  	d$other<-d$other+d$V2
  	d<-d %>% select(-V2)
  	print(j)
	}

	d$out<-d$other/(d$count+d$other)

	d<-d %>% filter(out<.01)
	#We're saving everything with less than 1% out of type overlaps, we'll filter more strictly in the next step
	write.table(d,sv,quote=F,row.names = F,sep="\t")
	print("YES!!!!!!")
	print(nrow(d))
} else {
	print("No need!")
}
