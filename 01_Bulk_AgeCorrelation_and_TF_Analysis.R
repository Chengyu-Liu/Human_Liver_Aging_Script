library(tidyverse)
library(dplyr)
library(reshape)
library(reshape2)
library(data.table)

library(edgeR)
library(DESeq2)
library(limma)

library(Seurat)
library(monocle)
library(SCENIC)
library(RcisTarget)

library(future)

use_condaenv("tensorflow_env", conda = "/opt/homebrew/anaconda3/bin/conda", required = TRUE)
plan("multicore", workers = 12)
options(future.globals.maxSize = 40*1024^3)


bulk_dir <- '/path_of_mapping_data_rpkm_dictionary/'
data <- read.csv(paste0(bulk_dir,'HealthySample_RNAseq_rpkm_normalization_scaled.csv'),row.names = 1)

# Input data description:
# `data` is a matrix of normalized and scaled RPKM values derived from bulk RNA-seq
# of healthy human liver samples. Rows correspond to individual samples, columns
# correspond to genes, and matrix entries represent gene expression levels.
# `ages` is a numeric vector containing the chronological ages of the corresponding
# samples, ordered identically to the rows of `data`.

# 01 Gene-age_correlation in healthy sample--------------------------------
ages <- sapply(rownames(data),function(vec){as.numeric(strsplit(vec,'_')[[1]][2])}  )
identical(names(ages),rownames(data))
correlations <- sapply(data, function(feature) cor(feature, ages))

features <- data
features$age <- ages
data <- features[,c(ncol(features),1:(ncol(features)-1))]
cor_results <- apply(data[, -1], 2, function(x) cor.test(data$age, x))
p_values <- sapply(cor_results, function(x) x$p.value)
head(p_values)
p_values <- p_values[names(correlations)]

df_cor <- data.frame(gene=names(correlations),correlation=correlations,pval=p_values)
df_cor$drc <- ifelse(df_cor$correlation> 0&df_cor$pval<0.05,'U',
                         ifelse(df_cor$correlation< 0&df_cor$pval<0.05,'D','N'))

#02 bulk transcritptional factors---------------------
# `df_cor` contains the results of the age–gene correlation analysis computed
df_cor <- df_cor[df_cor$pval<0.05,]
rownames(df_cor) <- df_cor$gene
df_cor_U <- df_cor[df_cor$correlation>0&df_cor$pval<0.05,]
df_cor_D <- df_cor[df_cor$correlation<0&df_cor$pval<0.05,]

srt_obj <- CreateSeuratObject(data)

for (drc in c('U','D')) {
  try({
    deg_drc <- get(paste0('df_cor_',drc))
    exp.mat <- LayerData(srt_obj,assay="RNA",layer="counts")
    
    genes <- intersect(unique(as.character(deg_drc$gene)),rownames(exp.mat))
    exp.mat <- as.matrix(exp.mat)
    exp.mat <- as.matrix(exp.mat[genes,])
    print(paste0('----------------------',drc,'_data loaded','----------------------'))
    
    TF_drc_dir=paste0(TF_dir,'bulk_',drc,'/')
    dir.create(TF_drc_dir)
    setwd(TF_drc_dir)
    org <- 'hgnc' 
    data(list='motifAnnotations_hgnc_v9',package='RcisTarget')
    motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
    dbDir <- 'data_base/TF_feature/hg19' # RcisTarget databases location
    myDatasetTitle <- paste0(drc,' DEGs') # choose a name for your analysis
    data(defaultDbNames)
    #data(list='motifAnnotations_hgnc_v9',package='RcisTarget')
    motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=1)  ##nCores为10??????十???叱?
    saveRDS(scenicOptions, file='int/scenicOptions.Rds')
    print(paste0('----------------------',drc,'_data initialized','----------------------'))
    
    # Gene filter
    genesKept <- geneFiltering(exp.mat, scenicOptions=scenicOptions,
                               minCountsPerGene=3*.01*ncol(exp.mat),
                               minSamples=ncol(exp.mat)*.01)
    exprMat_filtered <- exp.mat[genesKept, ]
    dim(exprMat_filtered)
    print(paste0('----------------------',drc,'_gene filtered','----------------------'))
    
    # Run Genie3
    runCorrelation(exprMat_filtered, scenicOptions)
    #judgement-------------------------------------------------------#
    judgement = tryCatch({
      #correct pipeline
      runGenie3(exprMat_filtered, scenicOptions, nParts = 10)
      print(paste0('----------------------',drc,'_Genie3 finished','----------------------'))
    },  error=function(e) e
    )
    if(inherits(judgement, 'error')) {
      print(paste0(drc,'_error'))
      next}
    #---------------------------------------------------------------#
    # Run the remaining
    scenicOptions@settings$verbose <- TRUE
    scenicOptions@settings$nCores <- 1
    scenicOptions@settings$seed <- 5
    #judgement-------------------------------------------------------#
    judgement = tryCatch({
      #correct pipeline
      runSCENIC_1_coexNetwork2modules(scenicOptions)
      runSCENIC_2_createRegulons(scenicOptions)
      runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
      print(paste0('----------------------',drc,'_',i,'_data completed','----------------------'))
    },  error=function(e) e
    )
    if(inherits(judgement, 'error')) {
      print(paste0(drc,'_',i,'_error'))
      next}
    #---------------------------------------------------------------#
  })
  
}




