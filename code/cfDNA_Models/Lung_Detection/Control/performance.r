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

dat %>% group_by(model) %>% summarize(r=round(roc(type,score)$auc[1],2))


write.csv(dat,"../Results/Cross_Validation_scores.csv")

