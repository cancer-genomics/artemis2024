library(data.table)
library(tidyverse)
library(devtools)
load_all("useful.stuff.aa")           # Load



features<-list.files("./counts/features",full.names=T)
dat<-squish(features)

colnames(dat)<-c("feature","count","id")

counts<-list.files("./read_counts",full.names=T)
counts<-rbindlist(lapply(counts, fread,header=F))
counts$V2<-rep(c("id","norm"),nrow(counts)/2)
counts<-cbind(counts %>% filter(V2=="id"),counts %>% filter(V2=="norm"))
colnames(counts)<-c("id","spare","norm","spare2") 
counts<-counts%>% select(id,norm)
dat<-inner_join(dat,counts,by="id")
dat$norm<-as.numeric(dat$norm)
dat$n_norm<-dat$count/(dat$norm/1000000)
dat<-dat[,fam.scaled:=scale(n_norm), by=id]
dat<-dat %>% mutate(set=sapply(str_split(feature,"#"),"[",2))
dat<-dat %>% mutate(set=if_else(is.na(set),"Centromere",set))

ref<-dat %>% group_by(id,set) %>% summarize(group_n=sum(count))
d<-inner_join(dat,ref,by=c("id","set"))
d$prop<-d$count/d$group_n
d<-d[,prop.scaled:=scale(prop), by=id]

data<-d %>% select(feature,n_norm,id) %>% spread(key=feature,value=n_norm)
write.csv(data,"artemis.csv")









