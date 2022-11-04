suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(pheatmap) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(plyr)
  library(dplyr)
  library(scuttle)
  library(Seurat)
  library(zellkonverter)
  library(SingleCellExperiment)
})
# take pipecomp results file path as input
args = commandArgs(trailingOnly = TRUE)
cluster.path <- args[[1]]
downsampleType <- args[[3]]

split_cluster <- function(x, scePath){
  print('clusters path')
  print(x)
  clusters <- readRDS(x)
  
  name <- strsplit(x, "pipecomp_outputs/")[[1]][[2]]
  name <- strsplit(name, "/")[[1]][[1]]
  print('name')
  print(name)
  print('scePath')
  print(scePath)
  sce <- readRDS(scePath)
  list.names <- c()
  for (i in 1:length(clusters)){
    list.names[i] <- names(clusters)[i]
    print("run name")
    run.name <- gsub(";", "+", names(clusters)[i])
    print(run.name)
    #file <- paste0(paste0(name, "/"), i, ".csv")
    file <- paste("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/clusters", 
                  name, 
                  downsampleType,
                  paste0(i, ".csv"),
                  sep="/")
    h5ad_filename <- paste0(i, ".h5ad")
    h5ad_filename <- paste("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/H5AD_files", 
                           name, 
                           downsampleType,
                           h5ad_filename,
                           sep="/")
    print('h5ad output name')
    if (is.na(clusters[i])){
      write.csv(c(run.name, "no clusters"), file)
      write.csv(c(run.name, "no clusters"), h5ad_filename)
    }
    else{
      write.csv(clusters[i], file)
      filt_sce <- sce[,rownames(as.matrix(clusters[[i]]))]
      print("subset ok")
      filt_sce <- logNormCounts(filt_sce)
      
      print('seurat')
      seu <- as.Seurat(filt_sce, counts="counts", data="logcounts")
      var.genes <- FindVariableFeatures(seu)
      var.names <- head(VariableFeatures(var.genes), 500)
      print(var.names)
      filt_sce <- filt_sce[var.names,]

      zellkonverter::writeH5AD(filt_sce, h5ad_filename)
    }
    
  }
  
  
  
}

split_cluster(cluster.path, scePath = args[[2]])