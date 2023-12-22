
library(data.table)
library(tidyverse)
library(devtools)

dat<-fread("ucsc-t2t-repeat-masker.bed")

data<-dat %>% select(`#chrom`,chromStart,chromEnd,name,strand)

data$subfam<-sapply(str_split(data$name,"#"),"[",2)
data$fam<-sapply(str_split(data$subfam,"/"),"[",1)

#there are 16 fams, I am going to exclude Simple_repeat,Low_complexity and Unknown
#The Simplerepeats are short 2-7mers, and their identification relies on many many 
#repeats (1000s) which is hard in shortread. 
#Single occurrences by chance will be high.
#Unknown will be hard to know if it is good biology or technical artifact
#Low complexity is like 1-2mers in many repeats, not good with short read
data<-data %>% filter(!fam %in% c("Unknown","Simple_repeat","Low_complexity"))
library(GenomicRanges)
data$chrom<-data$`#chrom`

s<-makeGRangesFromDataFrame(data,keep.extra.columns=T)
l<-data.frame(disjoin(s))
l<-l %>% group_by(seqnames,strand)%>%summarize(n=sum(width))
l$type<-"full_genome"
l$group<-"full_genome"
size<-l

names<-unique(data$fam)
for (i in 1:length(names)) {
  s<-makeGRangesFromDataFrame(data %>% filter(fam==names[i]),keep.extra.columns=T)
  l<-data.frame(disjoin(s)) %>% group_by(seqnames,strand)%>%summarize(n=sum(width))
  l$group<-names[i]
  l$type<-"fam"
  size<-rbind(size,l)
  print(i)

}

names<-unique(data$subfam)
for (i in 1:length(names)) {
  s<-makeGRangesFromDataFrame(data %>% filter(subfam==names[i]),keep.extra.columns=T)
  l<-data.frame(disjoin(s)) %>% group_by(seqnames,strand)%>%summarize(n=sum(width))
  l$group<-names[i]
  l$type<-"subfam"
  size<-rbind(size,l)
  print(i)
}

names<-unique(data$name)
for (i in 1:length(names)) {
  s<-makeGRangesFromDataFrame(data %>% filter(name==names[i]),keep.extra.columns=T)
  l<-data.frame(disjoin(s)) %>% group_by(seqnames,strand)%>%summarize(n=sum(width))
  l$group<-names[i]
  l$type<-"repeat_type"
  size<-rbind(size,l)
  print(i)
}

size$mb<-size$n/1000000
write.csv(size,"RepeatFamilySizes.csv")
