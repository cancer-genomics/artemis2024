library(plyr)
library(dplyr)
library(tidyverse)
library(caret)
library(recipes)
library(devtools)
library(data.table)
load_all("/dcs04/scharpf/data/annapragada/useful.stuff.aa")
library(pROC)
cohort_name<-"LUCAS"
#read in the cohort
train<-fread(paste0("../Cohorts/",cohort_name,"_Train.csv"),header=T) %>% select(-V1)
test<-fread(paste0("../Cohorts/",cohort_name,"_Test.csv"),header=T) %>% select(-V1)
#test<-fread(paste0("../Cohorts/",cohort_name,fold,"_Test.csv"),header=T)%>% select(-V1)
train_meta<-train
test_meta<-test

glmnetGrid <- expand.grid(
    alpha = 1,
    lambda = 10^seq(-5, -1, length.out = 100))

set.seed(1234)
ctrl_all <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 10,
                     verboseIter = TRUE,
                     savePredictions="final",
                     classProbs=TRUE,
                     index=createMultiFolds(train$type, 5, 10),
                     summaryFunction = twoClassSummary)

#######Train all the SSLs

#SSL 2 - Coverage GBM
model_name<-"Cov_GBM"
recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-starts_with("cov")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "gbm",
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-train_preds
preds_val<-test_preds
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

#######
#SSL 3 - Ratios alone
model_name<-"Ratios_ssl"
recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-starts_with("ratio"),-type) %>%
    step_pca(starts_with("ratio"), prefix = "ratio_pc_",threshold=.9)     %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

#######
#SSL 4 - zscores alone
model_name<-"zscores_ssl"
recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-starts_with("zscore")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

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

recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

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

recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

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

recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

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

recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

#SSL 9 - RNA_TE
model_name<-"RNA_TE_ssl"
e<-fread("/dcs04/scharpf/data/annapragada/mced_paper/Manuscript/code/modeling/Fixed_Versions_HiSeq/Expected.csv")
e<-e %>% filter(total_kmers>1000)

e$fam<-sapply(str_split(e$feature,"#"),"[",2)
e$fam<-sapply(str_split(e$fam,"_"),"[",1)
e<-e %>% mutate(fam=if_else(is.na(fam),"Satellite",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("rRNA","snRNA","scRNA","tRNA","srpRNA"),"RNA/DNA Elements",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("DNA","DNA?","RC","Retroposon"),"RNA/DNA Elements",fam))

e<-e %>% filter(fam=="RNA/DNA Elements")

recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-e$feature) 

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

#SSL 10 - Epi
model_name<-"Epi_ssl"
recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-type,-starts_with("states"),-starts_with("H3K"),-starts_with("H4K")) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())
colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-rbind(train_preds,preds_cv)
preds_val<-rbind(test_preds,preds_val)

##### Format the SSL features into matrix
preds_cv<-preds_cv %>% spread(key=model,value=score)
preds_val<-preds_val %>% spread(key=model,value=score)

train<-inner_join(train,preds_cv %>% select(-type),by="id")
test<-inner_join(test,preds_val %>% select(-type),by="id")
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

###Train Ensembles
#ARTEMIS_only
model_name<-"ARTEMIS_Ensemble"
recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-LINE_ssl,-Epi_ssl,-LTR_ssl,-RNA_TE_ssl,-id,-type,-Sat_ssl,-SINE_ssl)

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
test_preds<-get_test_preds(test,model1)
test_preds$model<-model_name
ensemble_res<-test_preds
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))

### save these scores for training other ensembles below
artemis_score_test<-test_preds
artemis_score_train<-get_cv_preds(train,model1)
artemis_score_train<-artemis_score_train %>% dplyr::rename("Artemis_Score"="score") %>% select(-type)
artemis_score_test<-artemis_score_test %>% dplyr::rename("Artemis_Score"="score")  %>% select(-type,-model)


#Set up the other kind of ensemble
train<-inner_join(train,artemis_score_train,by="id")
test<-inner_join(test,artemis_score_test,by="id")

#1 ARTEMIS score, trained with all the DELFI components
model_name<-"ARTEMIS_single_DELFI_SSLs_Ensemble"
recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-Artemis_Score,-id,-type,-Cov_GBM,-Ratios_ssl,-zscores_ssl)

colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
test_preds<-get_test_preds(test,model1)
test_preds$model<-model_name
ensemble_res<-rbind(test_preds,ensemble_res)
saveRDS(model1,paste0("../Locked_Models/",model_name,".rds"))



preds_val<-preds_val %>% gather(key=model,value=score,-id,-type)
preds_val$set<-"SSL"
ensemble_res$set<-"Ensemble"
res<-rbind(ensemble_res,preds_val)



write.csv(res,"../Results/Test_set_scores.csv")
















