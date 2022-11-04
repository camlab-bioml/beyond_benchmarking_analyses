suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(Matrix) # CRAN
  library(ggplot2)
  library(caret)
  library(glmnet)
  library(ggforce)
  library(SingleCellExperiment)
  library(pipeComp)
  library(scuttle)
})

args <- commandArgs(trailingOnly=TRUE)
sce_filename <- args[[1]]
sce <- readRDS(sce_filename)

source(system.file("extdata", "scrna_alternatives.R", package="pipeComp"))

# Downsample counts
counts.downsample <- downsampleMatrix(counts(sce), prop = 0.5)
counts.downsample <- Matrix(counts.downsample, sparse=TRUE)
sce.downsample <- SingleCellExperiment(list(counts = counts.downsample))

# computing QC metrics and metadata
sce.downsample <- addPerCellQC(sce.downsample)
sce.downsample <- addPerFeatureQC(sce.downsample)
sce.downsample <- add_meta(sce.downsample)
sce.downsample$phenoid <- sce$phenoid

# renaming some colData to match what's used in pipeComp
pct_counts_in_top_50_features <- sce.downsample$pct_counts_top_50_features
sce.downsample$pct_counts_in_top_50_features <- pct_counts_in_top_50_features
sce.downsample$pct_counts_top_50_features <- NULL
sce.downsample$pct_counts_Mt <- sce.downsample$pct_Mt
sce.downsample$pct_Mt <- NULL
sce.downsample$pct_counts_in_top_20_features <- sce.downsample$percent.top_20
sce.downsample$percent.top_20 <- NULL

# save to RDS
filename <- args[[2]]
saveRDS(sce.downsample, filename)

