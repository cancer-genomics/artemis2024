args <- commandArgs(trailingOnly = TRUE)

fold <- args[1]


library(plyr)
library(dplyr)
library(tidyverse)
library(caret)
library(recipes)
library(devtools)
library(data.table)
load_all("/dcs04/scharpf/data/annapragada/useful.stuff.aa")
library(pROC)
cohort_name<-"Cristiano"
#fold="fold10"
#read in the cohort
train<-fread(paste0("../Cohorts/",cohort_name,"_Train.csv"),header=T) %>% select(-V1)
fold=as.numeric(fold)
test<-train[fold]
train<-train %>% filter(!id %in% test$id)
#test<-fread(paste0("../Cohorts/",cohort_name,fold,"_Test.csv"),header=T)%>% select(-V1)
train_meta<-train
test_meta<-test

set.seed(1234)
ctrl_all <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 1,
                     verboseIter = TRUE,
                     savePredictions="final",
                     classProbs=TRUE,
                     index=createMultiFolds(train$Tumor_Type, 5, 1))

#######Train all the SSLs
#SSL 1 -- Mathios DELFI
model_name<-"Mathios_DELFI"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-starts_with("ratio"),-Tumor_Type,-starts_with("zscore"),-starts_with("cov")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)



train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-train_preds
preds_val<-test_preds
#######
#SSL 2 - Coverage GBM
model_name<-"Cov_GBM"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-starts_with("cov")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)

train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)

#######
#SSL 3 - Ratios alone
model_name<-"Ratios_ssl"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-starts_with("ratio"),-Tumor_Type) %>%
    step_nzv(all_predictors())

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)

#######
#SSL 4 - zscores alone
model_name<-"zscores_ssl"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-starts_with("zscore")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)

#######
#SSL 5 - LTR
model_name<-"LTR_ssl"
e<-fread("/dcs04/scharpf/data/annapragada/mced_paper/Manuscript/code/modeling/Fixed_Versions_HiSeq/Expected.csv")
e<-e %>% filter(total_kmers>1000)

e$fam<-sapply(str_split(e$feature,"#"),"[",2)
e$fam<-sapply(str_split(e$fam,"_"),"[",1)
e<-e %>% mutate(fam=if_else(is.na(fam),"Satellite",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("rRNA","snRNA","scRNA","tRNA","srpRNA"),"RNA/DNA Elements",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("DNA","DNA?","RC","Retroposon"),"RNA/DNA Elements",fam))

e<-e %>% filter(fam=="LTR")

recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)

#######
#SSL 6 - LINE
model_name<-"LINE_ssl"
e<-fread("/dcs04/scharpf/data/annapragada/mced_paper/Manuscript/code/modeling/Fixed_Versions_HiSeq/Expected.csv")
e<-e %>% filter(total_kmers>1000)

e$fam<-sapply(str_split(e$feature,"#"),"[",2)
e$fam<-sapply(str_split(e$fam,"_"),"[",1)
e<-e %>% mutate(fam=if_else(is.na(fam),"Satellite",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("rRNA","snRNA","scRNA","tRNA","srpRNA"),"RNA/DNA Elements",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("DNA","DNA?","RC","Retroposon"),"RNA/DNA Elements",fam))
e<-e %>% filter(fam=="LINE")

recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)


#######
#SSL 7 - Sat
model_name<-"Sat_ssl"
e<-fread("/dcs04/scharpf/data/annapragada/mced_paper/Manuscript/code/modeling/Fixed_Versions_HiSeq/Expected.csv")
e<-e %>% filter(total_kmers>1000)

e$fam<-sapply(str_split(e$feature,"#"),"[",2)
e$fam<-sapply(str_split(e$fam,"_"),"[",1)
e<-e %>% mutate(fam=if_else(is.na(fam),"Satellite",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("rRNA","snRNA","scRNA","tRNA","srpRNA"),"RNA/DNA Elements",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("DNA","DNA?","RC","Retroposon"),"RNA/DNA Elements",fam))
e<-e %>% filter(fam=="Satellite")

recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)

#SSL 8 - SINE
model_name<-"SINE_ssl"
e<-fread("/dcs04/scharpf/data/annapragada/mced_paper/Manuscript/code/modeling/Fixed_Versions_HiSeq/Expected.csv")
e<-e %>% filter(total_kmers>1000)

e$fam<-sapply(str_split(e$feature,"#"),"[",2)
e$fam<-sapply(str_split(e$fam,"_"),"[",1)
e<-e %>% mutate(fam=if_else(is.na(fam),"Satellite",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("rRNA","snRNA","scRNA","tRNA","srpRNA"),"RNA/DNA Elements",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("DNA","DNA?","RC","Retroposon"),"RNA/DNA Elements",fam))
e<-e %>% filter(fam=="SINE")

recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)

train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)

#SSL 9 - SINE
model_name<-"RNA_TE_ssl"
e<-fread("/dcs04/scharpf/data/annapragada/mced_paper/Manuscript/code/modeling/Fixed_Versions_HiSeq/Expected.csv")
e<-e %>% filter(total_kmers>1000)

e$fam<-sapply(str_split(e$feature,"#"),"[",2)
e$fam<-sapply(str_split(e$fam,"_"),"[",1)
e<-e %>% mutate(fam=if_else(is.na(fam),"Satellite",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("rRNA","snRNA","scRNA","tRNA","srpRNA"),"RNA/DNA Elements",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("DNA","DNA?","RC","Retroposon"),"RNA/DNA Elements",fam))

e<-e %>% filter(fam=="RNA/DNA Elements")

recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)

train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)


#SSL 10 - Epi
model_name<-"Epi_ssl"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-Tumor_Type,-starts_with("states"),-starts_with("H3K"),-starts_with("H4K")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())
colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)

train_preds<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
train_preds<-inner_join(train_preds,train_id,by="rowIndex")
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
train_preds<-train_preds %>% select(-rowIndex)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)


######
##### Format the SSL features into matrix
preds_cv<-preds_cv %>% gather(key=score_type,value=score,-pred,-id,-Tumor_Type,-model)
preds_cv$name<-paste0(preds_cv$model,"_",preds_cv$score_type)
preds_cv<-preds_cv %>% select(-score_type,-model,-pred)%>%spread(key=name,value=score)

preds_val_res<-preds_val
preds_val<-preds_val %>% gather(key=score_type,value=score,-pred,-id,-Tumor_Type,-model)
preds_val$name<-paste0(preds_val$model,"_",preds_val$score_type)
preds_val<-preds_val %>% select(-score_type,-model,-pred)%>%spread(key=name,value=score)


train<-inner_join(train,preds_cv %>% select(-Tumor_Type),by="id")
test<-inner_join(test,preds_val %>% select(-Tumor_Type),by="id")

###Train Ensembles
#ARTEMIS_only
model_name<-"ARTEMIS_Ensemble"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-starts_with("LINE_ssl"),-starts_with("Epi_ssl"),-starts_with("LTR_ssl"),
      -starts_with("RNA_TE_ssl"),-id,-Tumor_Type,-starts_with("Sat_ssl"),-starts_with("SINE_ssl"))

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
test_preds$model<-model_name
ensemble_res<-test_preds

### save these scores for training other ensembles below
artemis_score_test<-test_preds

artemis_score_train<-model1$pred %>% select(rowIndex,pred,Breast,Cholangiocarcinoma,Colorectal,Gastric,Lung,Ovarian,Pancreatic)
train_id<-train %>% select(id,Tumor_Type) %>% mutate(rowIndex=1:n())
artemis_score_train<-inner_join(artemis_score_train,train_id,by="rowIndex")
artemis_score_train <- artemis_score_train %>% select(-rowIndex,-pred,-Tumor_Type)
colnames(artemis_score_train)<-paste0("Artemis_Score_",colnames(artemis_score_train))
artemis_score_train<-artemis_score_train %>% dplyr::rename("id"="Artemis_Score_id") 

artemis_score_test <- artemis_score_test %>% select(-model,-pred,-Tumor_Type)

colnames(artemis_score_test)<-paste0("Artemis_Score_",colnames(artemis_score_test))
artemis_score_test<-artemis_score_test %>% dplyr::rename("id"="Artemis_Score_id") 

#ARTEMIS+Mathios
model_name<-"ARTEMIS_Mathios_Ensemble"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-starts_with("LINE_ssl"),-starts_with("Epi_ssl"),-starts_with("LTR_ssl"),
      -starts_with("RNA_TE_ssl"),-id,-Tumor_Type,-starts_with("Sat_ssl"),-starts_with("SINE_ssl"),-starts_with("Mathios_DELFI"))

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
test_preds$model<-model_name
ensemble_res<-rbind(test_preds,ensemble_res)

#ARTEMIS+DELFI_ssls
model_name<-"ARTEMIS_DELFI_SSls_Ensemble"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-starts_with("LINE_ssl"),-starts_with("Epi_ssl"),-starts_with("LTR_ssl"),
      -starts_with("RNA_TE_ssl"),-id,-Tumor_Type,-starts_with("Sat_ssl"),-starts_with("SINE_ssl"),-starts_with("Cov_GBM"),
      -starts_with("Ratios_ssl"),-starts_with("zscores_ssl"))

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
test_preds$model<-model_name
ensemble_res<-rbind(test_preds,ensemble_res)

#ARTEMIS+raw ratios and zscores
model_name<-"ARTEMIS_DELFI_raw_Ensemble"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-starts_with("LINE_ssl"),-starts_with("Epi_ssl"),-starts_with("LTR_ssl"),
      -starts_with("RNA_TE_ssl"),-id,-Tumor_Type,-starts_with("Sat_ssl"),-starts_with("SINE_ssl"),-starts_with("ratio"),
      -starts_with("zscore"),starts_with("zscores_ssl"),starts_with("Ratios_ssl")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())


colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
test_preds$model<-model_name
ensemble_res<-rbind(test_preds,ensemble_res)


#Set up the other kind of ensemble
train<-inner_join(train,artemis_score_train,by="id")
test<-inner_join(test,artemis_score_test,by="id")


#1 ARTEMIS score_single, trained with Mathios DELFI
model_name<-"ARTEMIS_single_Mathios_Ensemble"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-starts_with("Artemis_Score"),-id,-Tumor_Type,-starts_with("Mathios_DELFI"))

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
test_preds$model<-model_name
ensemble_res<-rbind(test_preds,ensemble_res)


#1 ARTEMIS score, trained with all the DELFI components
model_name<-"ARTEMIS_single_DELFI_SSLs_Ensemble"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-starts_with("Artemis_Score"),-id,-Tumor_Type,-starts_with("Ratios_ssl"),-starts_with("Cov_GBM"),-starts_with("zscores_ssl"))

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
test_preds$model<-model_name
ensemble_res<-rbind(test_preds,ensemble_res)


#1 ARTEMIS score, trained with raw delfi features
model_name<-"ARTEMIS_single_DELFI_raw_Ensemble"
recipe_seq <- recipe(Tumor_Type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-starts_with("Artemis_Score"),-id,-Tumor_Type,-starts_with("ratio"),-starts_with("zscore"),
      starts_with("zscores_ssl"),starts_with("Ratios_ssl")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())



colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
test_preds<-predict(model1,test,type="prob")
test_preds<-cbind(test %>% select(id,Tumor_Type),test_preds,pred=predict(model1,test))
test_preds$model<-model_name
ensemble_res<-rbind(test_preds,ensemble_res)




preds_val_res$set<-"SSL"
ensemble_res$set<-"Ensemble"
res<-rbind(ensemble_res,preds_val_res)



write.csv(res,paste0("../CV_Scores/","CV_scores_",fold,".csv"))

