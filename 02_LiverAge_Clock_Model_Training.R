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
#-------------------------------------#
bulk_dir <- '/path_of_mapping_data_rpkm_dictionary/'
##01.1 Load data-------------------
data <- read.csv(paste0(bulk_dir,'HealthySample_RNAseq_rpkm_normalization_scaled.csv'),row.names = 1)

# Input data description:
# `data` is a matrix of normalized and scaled RPKM values derived from bulk RNA-seq
# of healthy human liver samples. Rows correspond to individual samples, columns
# correspond to genes, and matrix entries represent gene expression levels.
# `ages` is a numeric vector containing the chronological ages of the corresponding
# samples, ordered identically to the rows of `data`.


##01.2 Split data-------------------
split_data <- function(data, prop = 0.75) {
  set.seed(123)  # for reproducibility
  split <- initial_split(data, prop = prop)
  list(
    train_data = training(split),
    test_data = testing(split)
  )
}


##01.3 Update data to remove botton weight genes-------------------
update_training_data <- function(train_data,  gene_bottom_ls) {
  gene_to_remove <- intersect(colnames(train_data), gene_bottom_ls)
  train_data <- train_data[, !colnames(train_data) %in% gene_to_remove]
  list(train_data = train_data)
}

##01.4 Adjust the weight of genes-------------------
adjust_gene_high_weights <- function(model, train_data, 
                                     gene_high, increase_factor = 1.2) {
  weights <- model$get_weights()
  gene_indices <- match(gene_high, colnames(train_data))
  gene_indices <- na.omit(gene_indices)  # Remove NAs
  if (length(gene_indices) > 0) {
    weights[[1]][gene_indices, ] <- weights[[1]][gene_indices, ] * increase_factor
    model$set_weights(weights)
  }
  return(model)
}



##01.5 Define the hyperparameter space-------------------
param_grid <- list(
  optimizer = c('rmsprop'), 
  learning_rate = c( 1e-4,3e-4,1e-3,3e-3,1e-2),  
  activation = c('relu'), 
  l1_reg = c(1e-4,3e-4,1e-3,3e-3,1e-2)
)

##01.6 Main training function-------------------
#Early stop to avoid overfitting#
callbacks <- callback_early_stopping(monitor = "val_loss", patience = 60, restore_best_weights = TRUE)

#Train and evaluate the model#
train_and_evaluate <- function(train_data, train_ages, 
                               layer_dims, 
                               param_grid,
                               round_number,
                               epochs) {
  optimizer_ls <- param_grid$optimizer
  learning_rate_ls <- param_grid$learning_rate
  activation_ls <- param_grid$activation
  l1_reg_ls <- param_grid$l1_reg
  #l2_reg_ls <- param_grid$l2_reg
  df_para_tmp <- NULL
  i=1
  for (optimizer in optimizer_ls) {
    for (learning_rate in learning_rate_ls) {
      for (activation in activation_ls) {
        for (l1_reg in l1_reg_ls) {
          #for (l2_reg in l2_reg_ls) {
          try({
            params <- list(optimizer = optimizer, 
                           learning_rate = learning_rate,  
                           activation = activation, 
                           l1_reg = l1_reg,
                           l2_reg = l1_reg)
            
            # Construct model
            model <- keras_model_sequential()
            model %>%
              layer_dense(units = layer_dims[1], activation = params$activation, input_shape = ncol(train_data),
                          kernel_regularizer = regularizer_l1_l2(l1 = params$l1_reg, l2 = params$l2_reg)) %>%
              layer_dense(units = layer_dims[2], activation = params$activation,
                          kernel_regularizer = regularizer_l1_l2(l1 = params$l1_reg, l2 = params$l2_reg)) %>%
              layer_dense(units = layer_dims[3], activation = params$activation,
                          kernel_regularizer = regularizer_l1_l2(l1 = params$l1_reg, l2 = params$l2_reg)) %>%
              layer_dense(units = layer_dims[4], activation = params$activation,
                          kernel_regularizer = regularizer_l1_l2(l1 = params$l1_reg, l2 = params$l2_reg)) %>%
              layer_dense(units = 1)
            
            # Compile model
            optimizer_selected <- switch(params$optimizer,
                                         "adam" = optimizer_adam(learning_rate = params$learning_rate),
                                         "rmsprop" = optimizer_rmsprop(learning_rate = params$learning_rate),
                                         "sgd" = optimizer_sgd(learning_rate = params$learning_rate),
                                         "adamax" = optimizer_adamax(learning_rate = params$learning_rate))
            
            model %>% compile(
              loss = 'mean_absolute_error',
              optimizer = optimizer_selected,
              metrics = list('mean_absolute_error')
            )
            
            ##01.6.1 Adjust weights for high correlation genes & gene_high-------------------------------
            model <- adjust_gene_high_weights(model, train_data, gene_high)
            print(paste0('#----- Round',round_number,'-No.',i,' ---- Adjust high frequency top gene weights completed'))
            
            
            
            ##01.6.2 Train data-----------------------------
            history <- model %>% fit(
              train_data, train_ages,
              epochs = epochs,
              callbacks = list(callbacks),
              verbose = 0
            )
            print(paste0('#----- Round',round_number,'-No.',i,' ---- Training data completed'))
            
            train_loss <- round(history$metrics$mean_absolute_error[length(history$metrics$mean_absolute_error)],1)
            
            print(paste0('#----- Round',round_number,'-No.',i,' ---- train_MAE: ',train_loss))
            
            ##01.6.3 Extract and normalize weights----------------------
            first_layer_weights <- model$get_weights()[[1]]
            first_layer_mean_contributions <- rowMeans(abs(first_layer_weights), na.rm = TRUE)
            min_val <- min(first_layer_mean_contributions)
            max_val <- max(first_layer_mean_contributions) - min_val
            normalized_contributions <- (first_layer_mean_contributions - min_val) / max_val
            contributions <- setNames(as.list(normalized_contributions), colnames(train_data))
            
            ##01.6.4 Record training iteration---------------
            current_row <- data.frame(
              node_layers = toString(layer_dims),
              gene_count = ncol(train_data),
              optimizer = params$optimizer,
              learning_rate = params$learning_rate,
              activation = params$activation,
              l1_reg = l1_reg ,
              l2_reg = l1_reg ,
              train_MAE = train_loss,
              round = round_number,
              model=paste0('Model_',i)
            )
            ##01.6.5  merge gene weights----------------
            current_row <- cbind(current_row, contributions)
            df_para_tmp <- rbind(df_para_tmp, current_row)
            if (train_loss<=5) {
              model %>% save_model_hdf5(paste0(model_wd,'Model_Round',round_number,'-No.',i,
                                               '_trainMAE-',train_loss,
                                               '_Gn',n_genes,
                                               '.hdf5'))
              saveRDS(train_data,paste0(model_wd,'Model_Round',round_number,'-No.',i,
                                        '_trainMAE-',train_loss,
                                        '_Gn',n_genes,
                                        '_trainData.rds'))
            }
            
            
            i=i+1
          })
          
          #}
        }
      }
    }
  }
  
  return(df_para_tmp)
}

##01.7 Main training cycle-------------------
# data_frames <- initialize_data_frames()
split <- split_data(data)
train_data <- as.matrix(split$train_data)
train_ages <- ages[rownames(train_data)]
layer_dims <- c(512, 512, 256, 256)
epochs <- 500

n_genes <- ncol(train_data)
gene_top_tt <- NULL
gene_bottom_ls <- NULL
gene_high <- NULL

df_para <- NULL
round_number <- 1
found_satisfactory_model <- FALSE

while (found_satisfactory_model == F | n_genes<10) {
  print(paste0('#----- Round',round_number,' initializated ------ GeneNumber: ',n_genes,' -----#'))
  df_para_tmp <- train_and_evaluate(train_data, train_ages, 
                                    layer_dims, 
                                    param_grid,
                                    round_number,
                                    epochs)
  df_para_tmp <- df_para_tmp[is.na(df_para_tmp$train_MAE)==F ,]
  
  
  df_para <- rbind(df_para,df_para_tmp[,c(1:11)])
  write.csv(df_para,paste0(model_wd,'df_para.csv'),row.names = F)
  
  ###01.7.1 Selecting top 5 best models--------------------
  sorted_results <- df_para_tmp %>% arrange(train_MAE)
  num <- nrow(sorted_results[sorted_results$train_MAE<=5,])
  if (num<=5) {
    best_models <- head(sorted_results, 5 )
  }else {
    best_models <- sorted_results[sorted_results$train_MAE<=5,]
  }
  
  
  ##01.7.2 Calculate the average weight of the selected models--------------------
  gene_contributions <- best_models[,c(12:ncol(best_models))]
  avg_contributions <- colMeans(gene_contributions)
  
  
  ##01.7.3 Sequence the gene contribution-------------------
  gene_order <- order(avg_contributions, decreasing = T)
  ordered_genes <- data.frame(gene=names(avg_contributions)[gene_order],order=1:length(gene_order))
  
  df_para_epoh_tmp <- best_models %>% 
    summarise(
      node_layers = first(node_layers), 
      gene_count = first(gene_count), 
      train_MAE_mean = mean(train_MAE)
    ) %>%
    mutate(round = round_number)
  
  df_para_epoh_tmp <- cbind(df_para_epoh_tmp, ordered_genes)
  
  ##01.7.4 Test model quality--------------
  try({
    df_catch <- sorted_results[sorted_results$train_MAE<=5,]
  })
  if (nrow(df_catch)>0) {
    cat("Found a satisfactory model with train MAE: ", df_catch$train_MAE, "\n")
    found_satisfactory_model <- TRUE
    break  
  }
  
  ##01.7.5 Select the top and bottom genes--------------
  gene_top_ls <- df_para_epoh_tmp[df_para_epoh_tmp$order<=max(df_para_epoh_tmp$order)*0.05,]$gene
  gene_bottom_ls <- df_para_epoh_tmp[df_para_epoh_tmp$order>=max(df_para_epoh_tmp$order)*0.95,]$gene
  
  
  gene_top_tt <- c(gene_top_tt,gene_top_ls)
  gene_high <- names(table(gene_top_tt)[table(gene_top_tt)>=5])
  
  
  ##01.7.6 Updata train-----------------
  updated_data <- update_training_data(train_data,  gene_bottom_ls)
  train_data <- updated_data$train_data
  
  n_genes <- ncol(train_data)
  round_number <- round_number+1
  save.image(paste0(model_wd,'model_tmp.Rdata'))
  
}

