suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(SingleCellExperiment)
  library(pcaMethods)
  library(stats)
})

args <- commandArgs(trailingOnly = TRUE)
trainTestSplit <- readRDS(args[[length(args)]])

# Get average expression for each dataset
means <- args[-length(args)] %>%
  map(readRDS) %>%
  map(counts) %>%
  map(rowMeans) %>%
  map(as.data.frame)%>%
  map(rownames_to_column,var="gene") %>%
  purrr::reduce(merge, by="gene", all=TRUE)%>%
  column_to_rownames(var="gene")

# Get experiment names 
expNames <- strsplit(args[-length(args)], split="/") %>%
  map(grep, value=TRUE, pattern="^E") %>%
  map(2) %>%
  map(strsplit, split=".RDS") %>%
  map(1)

print(expNames)
colnames(means) <- expNames
#print(means)

# Split into train/test
trainMeans <- means %>%
  select(starts_with(trainTestSplit$train))

testMeans <- means %>% 
  select(starts_with(trainTestSplit$test))

trainMeans <- t(trainMeans)
testMeans <- t(testMeans)

# print(trainMeans)
# print(testMeans)

print(dim(trainMeans))
print(dim(testMeans))

trainMeans[is.na(trainMeans)] <- 0
testMeans[is.na(testMeans)] <- 0

print(any(is.na(trainMeans)))
print(any(is.na(testMeans)))

pcaTrain <- pca(as.matrix(trainMeans), nPcs=20, seed=123)
print(pcaTrain)
pcaTest <- predict(pcaTrain, newdata = as.matrix(testMeans))

trainExprScores <- pcaTrain@scores
testExprScores <- pcaTest$scores

saveRDS(pcaTrain, "/home/campbell/cfang/automl_scrna/results/gene_expr_pca_model_downsampled.RDS")
saveRDS(list(trainExpr=trainExprScores, testExpr=testExprScores), "/home/campbell/cfang/automl_scrna/data/avg_expr_pca_downsampled.RDS")