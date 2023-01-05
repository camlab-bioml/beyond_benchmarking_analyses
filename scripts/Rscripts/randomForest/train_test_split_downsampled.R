suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(stats)
  library(stringr)
})
args <- commandArgs(trailingOnly=TRUE)
output <- args[[1]]
designMat <- readRDS(args[[2]])
trainTestSplit <- readRDS(args[[3]])

print(unique(designMat$name))
print(colnames(designMat))
trainNames <- trainTestSplit$train
testNames <- trainTestSplit$test

designMatTrain <- designMat %>%
  filter(apply(sapply(trainNames, FUN = grepl, .$name), MARGIN=1, FUN=any))

designMatTest <- designMat %>%
  filter(apply(sapply(testNames, FUN = grepl, .$name), MARGIN=1, FUN=any))

print(unique(designMatTest$name))
print(unique(designMatTrain$name))
print(intersect(designMatTrain$name, designMatTest$name))

split <- list(test=unique(designMatTest$name), train=unique(designMatTrain$name))
saveRDS(split, "/home/campbell/cfang/automl_scrna/data/train_test_split_downsampled.RDS")