args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])


library(data.table)
library(tidyverse)
library(devtools)
load_all("useful.stuff.aa")           # Load




f<-list.files("./denovo_kmer_type",pattern=".txt",full.names=T)

d<-fread(f[i])
sv<-gsub("./denovo_kmer_type","./kmers_non_masked",f[i])
sv<-gsub(".txt",".csv",sv)

if (!file.exists(sv)) {
	print("Running the code here!")
	test_dat<-fread("masked.txt")
	d<- d %>% rename("kmer"="V1")
	d<- d %>% rename("count"="V2")
	d<-d%>%filter(!kmer %in% test_dat$V1)
	write.csv(d,sv)
	print(i)
	print(sv)
}


