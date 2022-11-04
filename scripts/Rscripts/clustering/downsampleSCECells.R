suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(Matrix) # CRAN
  library(ggplot2)
  library(caret)
  library(glmnet)
  library(ggforce)
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly=TRUE)
sce_filename <- args[[1]]
sce <- readRDS(sce_filename)

nSamp50 <-floor(0.5*ncol(sce))
nSamp75 <- floor(0.75*ncol(sce))

set.seed(0)
sce50 <- sce[,sample(ncol(sce), size = nSamp50)]
sce75 <- sce[,sample(ncol(sce), size = nSamp75)]

dim(sce)
dim(sce50)
dim(sce75)

sce50_filename <- args[[2]]
sce75_filename <- args[[3]]
saveRDS(sce50, sce50_filename)
saveRDS(sce75, sce75_filename)

