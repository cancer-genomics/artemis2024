library(data.table)
library(tidyverse)
library(devtools)
load_all("useful.stuff.aa")           # Load

dat<-fread("ucsc-t2t-repeat-masker.bed")

data<-dat %>% select(`#chrom`,chromStart,chromEnd,name)

data$subfam<-sapply(str_split(data$name,"#"),"[",2)
data$fam<-sapply(str_split(data$subfam,"/"),"[",1)

#there are 16 fams, I am going to exclude Simple_repeat,Low_complexity and Unknown
#The Simplerepeats are short 2-7mers, and their identification relies on many many 
#repeats (1000s) which is hard in shortread. 
#Single occurrences by chance will be high.
#Unknown will be hard to know if it is good biology or technical artifact
#Low complexity is like 1-2mers in many repeats, not good with short read
data<-data %>% filter(!fam %in% c("Unknown","Simple_repeat","Low_complexity"))

stats<-data
stats$width<-stats$chromEnd-stats$chromStart
stats<-stats %>% group_by(fam) %>% summarize(total=sum(width))
stats$total_mb<-stats$total/1000000
stats<-stats %>% select(fam,total_mb)
#This leaves 57 subfams from 13 families, comprising 1287 types
#937 of these appear on the y chromosome and 1179 of them appear on the x chromosome

d1<-data %>% group_by(name,`#chrom`)%>%summarize(c=n())
d2<-data %>% group_by(name)%>%summarize(t=n())

d1<-inner_join(d1,d2,by="name")
d1$p<-d1$c/d1$t
d1 %>% filter(`#chrom`=="chrX" & p>.1) ### 106 appear on the x chr more than 10% of the time
d1 %>% filter(`#chrom`=="chrY" & p>.05) ### 77 appear on the y chr more than 10% of the time
#I'll still process these, but I should remember them!
#I could later remove kmers that are mapping a large proportion to chrX


names<-unique(data$name)
dir.create("bed_type")
for (i in 1:length(names)) {
  print(i)
  print(names[i])
  d<-data %>% filter(name==names[i])
  nm<-gsub("/","_",names[i])
  print(nm)
  write.table(d,paste0("./bed_type/",nm,".bed"),quote=F,row.names = F,sep="\t")
}