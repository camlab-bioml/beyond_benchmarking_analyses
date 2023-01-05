suppressPackageStartupMessages({
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(pipeComp) # Github
  library(purrr)
  library(scuttle)
  library(Seurat)
})
args <- commandArgs(trailingOnly = TRUE)
exps <- args[-c(1,2)]
downsampleType <- args[[2]]
outputFile <- args[[1]]

allColdata <- c()
expNames <- c()

for (i in 1:length(exps)){
  exp <- exps[[i]]
  print(exp)
  exp.name <- strsplit(exp, split="/")[[1]][[4]]
  print(exp.name)
  expNames <- c(expNames, exp.name)
  sce <- readRDS(exp)
  cd <- colData(sce)
  cd <- as.matrix(cd)
  cd <- apply(cd, 2, as.numeric)
  
  dataNames <- colnames(cd)
  dataNames <- dataNames[-which(dataNames %in% c("Assay", "phenoid"))]
  
  cd <- colMedians(cd[,-which(colnames(cd) %in% c("Assay", "phenoid"))])
  names(cd) <- dataNames
  
  cd <- c(cd, ncells=dim(sce)[[2]])
  cd <- c(cd, ngenes=dim(sce)[[1]])
  print(cd)
  allColdata[[i]] <- data.frame(cd)
  names(allColdata)[[i]] <- exp.name
}
print(allColdata)
allColdata <- do.call(cbind, allColdata)
colnames(allColdata) <- paste(expNames, downsampleType, sep="-")
out.name <- strsplit(exp, split=exp.name)[[1]][[1]]
allColdata <- t(allColdata)
allColdata <- as.data.frame(cbind(X=expNames, allColdata))
allColdata$downsampleType <- downsampleType
print(allColdata)
saveRDS(allColdata, outputFile)
