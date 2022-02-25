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

split_cluster <- function(x){
  clusters <- readRDS(x)
  
  name <- strsplit(x, "pipecomp_outputs/")[[1]][[2]]
  name <- strsplit(name, "/")[[1]][[1]]
  sce.name <- paste("/home/campbell/cfang/automl_scrna/data/experiments/", name, sep="")
  sce.name <- paste(sce.name, paste0(name, "-SCE.RDS"), sep="/")
  print(sce.name)
  sce <- readRDS(sce.name)
  
  list.names <- c()
  for (i in 1:length(clusters)){
    list.names[i] <- names(clusters)[i]
    run.name <- gsub(";", "+", names(clusters)[i])
    file <- paste0(paste0(name, "/"), i, ".csv")
    file <- paste0("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/clusters/", file)
    write.csv(clusters[i], file)
    
    filt_sce <- sce[,rownames(as.matrix(clusters[[i]]))]
    filt_sce <- logNormCounts(filt_sce)
    
    seu <- as.Seurat(filt_sce, counts="counts", data="logcounts")
    var.genes <- FindVariableFeatures(seu)
    var.names <- head(VariableFeatures(var.genes), 500)
    filt_sce <- filt_sce[var.names,]
    
    h5ad_filename <- paste0(paste0(name, "/"), i, ".h5ad")
    h5ad_filename <- paste0("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/H5AD_files/", h5ad_filename)
    print(h5ad_filename)
    zellkonverter::writeH5AD(filt_sce, h5ad_filename)
  }

  
  
}

split_cluster(cluster.path)