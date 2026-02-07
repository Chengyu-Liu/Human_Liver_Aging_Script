library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(patchwork)
library(ggridges)
library(SoupX)
library(cowplot)
library(reshape)
library(reshape2)
library(ggpubr)
library(ComplexHeatmap)
library(ggrepel)
library(scCustomize)
library(readr)
set.seed(5)
#----------------------------------------------------------------------------------------------#
# `sc_sp_ls` and `sn_sp_ls` store the sample identifiers for single-cell RNA-seq 
# and single-nucleus RNA-seq datasets, respectively.
# `sc_age_ls` and `sn_age_ls` contain the chronological ages corresponding to the
# single-cell RNA-seq and single-nucleus RNA-seq samples, respectively.

sample_list <- c(sc_sp_ls,sn_sp_ls)
sample_number <- length(sample_list)
age_list <- c(sc_age_ls,sn_age_ls)
#01.1 scRNAseq/snRNAseq data qc-----------------------------------
for (i in seq(1:sample_number)) {
  sp <- sample_list[i]
  if (sp%in%sc_sp_ls) {
    tmp <-  Read10X(paste0(mapping_dir[i],sp,'/outs/'))
    tmp <- CreateSeuratObject(counts = tmp, project = sp, min.cells = 3, min.features = 200)
    tmp$sample <- name_list[i]
    tmp$age<- age_list[i]
    tmp[["percent.MT"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
    tmp[["percent.RBL"]] <- PercentageFeatureSet(tmp, pattern = "^RP[SL]")
    assign(sample_list[i], tmp)
  }else{
    tmp <- Read_CellBender_h5_Mat(paste0(mapping_dir,sp,'/',sp,'_rawmatrix_filtered.h5'), use.names = F)
    tmp <- CreateSeuratObject(counts = tmp, project = sp, min.cells = 3, min.features = 200)
    tmp$sample <- name_list[i]
    tmp$age<- age_list[i]
    tmp[["percent.MT"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
    tmp[["percent.RBL"]] <- PercentageFeatureSet(tmp, pattern = "^RP[SL]")
    assign(sp, tmp)
  }
  
}
##01.2 set mito_percent---------------------------------
for (i in seq(1:sample_number)) {
  sp <- sample_list[i]
  if (sp%in%sc_sp_ls) {
    tmp <- get(sp)
    tmp <- subset(tmp, subset = percent.MT < 15 & nFeature_RNA > 600 & nFeature_RNA < 6000)
    tmp <- RenameCells(tmp, add.cell.id = name_list[i])
    assign(sample_list[i], tmp)
  }else {
    tmp <- get(sp)
    tmp <- subset(tmp, subset = percent.MT < 5 & nFeature_RNA > 200 & nFeature_RNA < 6000)
    tmp <- RenameCells(tmp, add.cell.id = name_list[i])
    assign(sp, tmp)
  }
  
}
##01.3 normalization --------------------------------
for (i in seq(1:sample_number)) {
  tmp <- get(sample_list[i])
  print(date())
  print(paste0(sample_list[i], ': SCTransform started'))
  tmp <- SCTransform(tmp, vars.to.regress = "percent.MT", verbose = FALSE,do.scale = T)
  print(date())
  print(paste0(sample_list[i], ': SCTransform finished'))
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_SCT.rds"))
  tmp <- RunPCA(tmp, verbose=F)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_bfPCR.rds"))
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
  tmp <- RunUMAP(tmp, dims = 1:20, verbose=F)
  tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:20)
  tmp <- FindClusters(tmp, res=res[i])
  tmp[["cluster"]] <- Idents(tmp)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_PCR.rds"))
  assign(sample_list[i], tmp)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}
##01.4 Identify doublets -----------------------------------
pK.df <- data.frame(matrix(nrow=0, ncol=2))
colnames(pK.df) <- c("Sample", "Optimal_pK")

for (i in c(1:sample_number)){
  sweep.res.list <- paramSweep(get(sample_list[i]), PCs = 1:20, sct = T, num.cores=16)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- arrange(bcmvn, desc(BCmetric))$pK[1]
  tmp <- data.frame(Sample=name_list[i], Optimal_pK=pK)
  pK.df <- rbind(pK.df, tmp)
  print(bcmvn)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}
write.csv(pK.df,file =paste0(qc_dir,ts,"_Optimal_pK.csv"), sep = "" )

sample_inforation <- data.frame(matrix(nrow=0, ncol=6))
colnames(sample_inforation) <- c("Sample", "Number","Doublet_prop","AvailableCellNumber","Cell_num_pre","Ratio_use")

for (i in seq(1:sample_number)) {
  tmp <- get(sample_list[i])
  pK.use <- as.numeric(as.character(pK.df$Optimal_pK[i]))
  homotypic.prop <- modelHomotypic(tmp@meta.data$cluster)
  
  cell_num_pre <- as.numeric(length(tmp@meta.data$orig.ident))
  ratio <- cell_num_pre/1000*0.008
  
  nExp_poi <- round(ratio*length(tmp@meta.data$orig.ident))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder(tmp, PCs = 1:20, pN = 0.25, pK = pK.use, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  tmp[["doublet"]] <- tmp[[paste("DF.classifications_0.25", pK.use, nExp_poi.adj, sep="_")]]
  prop <- nExp_poi.adj/length(tmp@meta.data$cluster)
  prop.tmp <- data.frame(Sample=name_list[i], Number=nExp_poi.adj, Doublet_prop=prop,
                         AvailableCellNumber=length(tmp@meta.data$doublet[tmp@meta.data$doublet=='Singlet']),
                         Cell_num_pre=cell_num_pre,Ratio_use=ratio)
  #
  sample_inforation <- rbind(sample_inforation, prop.tmp)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_doublets.rds"))
  assign(sample_list[i], tmp)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}

##01.5 filter out doublet -------------------------
for (i in seq(1:sample_number)){
  tmp <- get(sample_list[i])
  tmp <- subset(tmp, doublet=='Singlet')
  tmp <- SCTransform(tmp, verbose = FALSE, do.scale = T)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_final.rds"))
  print(paste0(i, ' completed', " (", i, "/" , sample_number ,")"))
}

##01.6 Integration ------------------------------
for (sp in sc_sp_ls) {
  tmp <- readRDS(file = paste0(qc_dir,sp,"_final.rds"))
  DefaultAssay(tmp) <- 'RNA'
  assign(sp,tmp)
}

int.list <- c()
for (i in sc_sp_ls) {
  tmp <- get(i)
  tmp2 <- list(tmp)
  int.list <- c(int.list,tmp2)
}

sc_srt <- merge(int.list[[1]], 
                    y = c(unlist(int.list[2:length(int.list)])),
                    project = "Human_liver")
DefaultAssay(sc_srt) <- 'RNA'
sc_srt <- NormalizeData(sc_srt,assay='RNA')

for (sp in sn_sp_ls) {
  tmp <- readRDS(file = paste0(qc_dir,sp,"_final.rds"))
  DefaultAssay(tmp) <- 'RNA'
  assign(sp,tmp)
}

int.list <- c()
for (i in sn_sp_ls) {
  tmp <- get(i)
  tmp2 <- list(tmp)
  int.list <- c(int.list,tmp2)
}

sn_srt <- merge(int.list[[1]], 
                y = c(unlist(int.list[2:length(int.list)])),
                project = "Human_liver")
DefaultAssay(sn_srt) <- 'RNA'
sn_srt <- NormalizeData(sn_srt,assay='RNA')


sc_srt <- FindVariableFeatures(sc_srt,assay='RNA')
sn_srt <- FindVariableFeatures(sn_srt,assay='RNA')
hvg_slt <- intersect(VariableFeatures(sc_srt),
                     VariableFeatures(sn_srt))
hvg_slt <- setdiff(hvg_slt,mt_genes)



int_sc <- SplitObject(sc_srt,split.by='sample')
int_sn <- SplitObject(sn_srt,split.by='sample')
int.list <- c(int_sc,int_sn)
merged_obj <- merge(int.list[[1]], 
                    y = c(unlist(int.list[2:length(int.list)])),
                    project = "Human_liver")
DefaultAssay(merged_obj) <- 'RNA'
merged_obj <- NormalizeData(merged_obj,assay='RNA')
merged_obj <- FindVariableFeatures(merged_obj,assay='RNA')
merged_obj <- ScaleData(merged_obj,assay='RNA',features=hvg_slt)
merged_obj_int <- RunPCA(merged_obj,features=hvg_slt)


mtd <- 'JointPCAIntegration'
merged_obj <- IntegrateLayers(
  object = merged_obj_int,
  method = get(mtd),
  orig.reduction = "pca", 
  new.reduction = "integrated_srt",
  verbose = F
)
merged_obj <- JoinLayers(merged_obj,assay='RNA')













