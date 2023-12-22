library(plyr)
library(dplyr)
library(tidyverse)
library(caret)
library(recipes)
library(devtools)
library(data.table)
load_all("/dcs04/scharpf/data/annapragada/useful.stuff.aa")
library(pROC)

f<-list.files("../CV_Scores",full.names=T)
dat<-squish(f)

dat %>% group_by(model,Tumor_Type) %>% dplyr::summarize(n=n(),c=sum(Tumor_Type==pred),p=c/n)


write.csv(dat %>% select(-V1),"../Results/Cross_Validation_scores.csv")

