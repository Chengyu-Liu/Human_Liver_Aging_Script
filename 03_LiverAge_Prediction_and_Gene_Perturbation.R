library(tidyverse)
library(dplyr)
library(reshape)
library(reshape2)
library(data.table)

library(reticulate)
library(tensorflow)
library(keras)
library(rhdf5)
library(caret)
library(rsample)
library(glmnet)
library(randomForest)

library(future)
reticulate::use_condaenv("R4.3.3")
Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")
Sys.setenv("TF_CPP_MIN_LOG_LEVEL" = "2")
options(future.globals.maxSize = 40*1024^3)

bulk_dir <- '/path_of_mapping_data_rpkm_dictionary/'
##01.1 Load data-------------------
data_all_sp <- read.csv(paste0(bulk_dir,'All_Sample_RNAseq_rpkm_normalization_scaled.csv'),row.names = 1)

# Input data description:
# `data_all_sp` is a matrix of normalized and scaled RPKM values derived from
# bulk RNA-seq of human liver samples, including both healthy and diseased individuals.
# Rows correspond to individual samples, columns correspond to genes, and matrix
# entries represent gene expression levels.
# `df_sample_info` is a data frame containing sample-level metadata, including
# individual identifiers, chronological age, and disease group annotation.


model <- load_model_hdf5('Model_healthy_liver_age_pred.hdf5')
df_sample_info_dse <- df_sample_info[df_sample_info$group=='disease',]
data_dse <- data_all_sp[df_sample_info_dse$sample,gene_clock]

df_sample_info_hty <- df_sample_info[df_sample_info$group=='healthy',]
data_hty <- data_all_sp[df_sample_info_hty$sample,gene_clock]

#02 Age prediction--------------------
##02.1 age prediction in healthy sample--------------
df_hty_pred <- prd_fn(data_hty,model)

##02.2 age prediction in disease sample--------------
df_dse_prid <- prd_fn(data_dse,model)

#03 Gene perturbation----------------
fit <- lm(formula = age_prd ~ age_ori, data = df_hty_pred)
df_hty_pred$dt_age <- df_hty_pred$age_prd-predict(fit,newdata = df_hty_pred)
df_hty_pred <- df_hty_pred[order(df_hty_pred$dt_age,decreasing = T),]


##03.1 acceleration-------
accelerator <- df_hty_pred[1:6,]$sample

df_pert_pro <- NULL
for (sp in accelerator) {
  data_slt <- data[c(sp,sp),colnames(train_data)]
  for (gn in colnames(data_slt)) {
    data_slt2 <- data_slt
    if (cor_df[cor_df$gene==gn,]$correlation>0) {
      data_slt2[1,gn] <- data_slt2[1,gn]-abs(data_slt2[1,gn]*0.25)
    }else {
      data_slt2[1,gn] <- data_slt2[1,gn]+abs(data_slt2[1,gn]*0.25)
    }
    new_predictions <- model %>% predict(data_slt2)
    dt_ag_tmp <- new_predictions[1,1]-new_predictions[2,1]
    df_tmp <- data.frame(gene=gn,pert=dt_ag_tmp,sample=sp)
    df_pert_pro <- rbind(df_pert_pro,df_tmp)
  }
  
}
df_pert_pro <- df_pert_pro%>%group_by(sample)%>%arrange(sample,desc(pert) )%>%mutate(ord=row_number())%>%as.data.frame()

df_pert_pro2 <- df_pert_pro%>%group_by(sample)%>%arrange(sample,pert )%>%mutate(ord=row_number())%>%as.data.frame()
df_pert_slt <- df_pert_pro2[df_pert_pro2$ord<=30,]
table(df_pert_slt$gene)[order(table(df_pert_slt$gene),decreasing = T)]


##03.2 deceleration--------
decelerator <- df_hty_pred[(nrow(df_hty_pred)-6):nrow(df_hty_pred),]$sample
df_pert_dly <- NULL
for (sp in decelerator) {
  data_slt <- data[c(sp,sp),colnames(train_data)]
  for (gn in colnames(data_slt)) {
    data_slt2 <- data_slt
    if (cor_df[cor_df$gene==gn,]$correlation>0) {
      data_slt2[1,gn] <- data_slt2[1,gn]+abs(data_slt2[1,gn]*0.25)
    }else {
      data_slt2[1,gn] <- data_slt2[1,gn]-abs(data_slt2[1,gn]*0.25)
    }
    new_predictions <- model %>% predict(data_slt2)
    dt_ag_tmp <- new_predictions[1,1]-new_predictions[2,1]
    df_tmp <- data.frame(gene=gn,pert=dt_ag_tmp,sample=sp)
    df_pert_dly <- rbind(df_pert_dly,df_tmp)
  }
  
}
df_pert_dly <- df_pert_dly%>%group_by(sample)%>%arrange(sample,pert)%>%mutate(ord=row_number())%>%as.data.frame()

df_pert_slt <- df_pert_dly[df_pert_dly$ord<=20,]
table(df_pert_slt$gene)[order(table(df_pert_slt$gene),decreasing = T)]