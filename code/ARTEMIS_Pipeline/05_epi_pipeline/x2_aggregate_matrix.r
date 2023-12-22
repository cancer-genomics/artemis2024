
library(tidyverse)
library(devtools)
load_all("useful.stuff.aa")
library(data.table)


f<-list.files("./epi_bins",full.names=T)

d<-squish(f)

setkey(d, id,size)
d[,ratio.cor:=short.cor/long.cor]
d[,ratio.scaled:=scale(ratio.cor), by=list(id,size)]
d[,cov.cor:=short.cor+long.cor]
d[,cov.scaled:=scale(cov.cor), by=list(id,size)]

epi_data<-d 
epi_data$bin<-paste0(epi_data$ref,"_",epi_data$bin)
epi_data<-epi_data%>%select(id,bin,cov.scaled) %>% spread(key=bin,value=cov.scaled)

write.csv(epi_data,"epi_bins.csv")

