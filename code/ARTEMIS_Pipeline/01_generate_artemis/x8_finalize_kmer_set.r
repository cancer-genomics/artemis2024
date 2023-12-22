args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

library(data.table)
library(tidyverse)
library(devtools)
load_all("useful.stuff.aa")           # Load

f<-list.files("./unique_kmer_type",pattern=".txt",full.names=T)
f<-f[i]
f2<-gsub("./unique_kmer_type","./kmers_non_masked",f)
f2<-gsub(".txt",".csv",f2)

sv<-gsub("./unique_kmer_type","./final_kmers",f)
sv<-gsub(".txt",".fasta",sv)

if (!file.exists(sv)) {
  print("Running the code here!")
  d<-fread(f)

  test_dat<-fread(f2)
  d<-d%>%filter(kmer %in% test_dat$kmer)
  d<-d %>% filter(out==0) #Here we ensure that they are unique to the repeat type (this table contains all with less than 1% out of family occurrences)
  d$type<-sapply(str_split(f,"/"),"[",3)
  d$type<-sapply(str_split(d$type,".txt"),"[",1)
  d$type<-paste0(">",d$type,"_",c(1:nrow(d)))
  d<-d %>% select(type,kmer)
  D <- do.call(rbind, lapply(seq(nrow(d)), function(i) t(d[i, ])))
  write.table(D, sv,row.names = FALSE, col.names = FALSE, quote = FALSE)

  print(i)
  print(sv)
}




