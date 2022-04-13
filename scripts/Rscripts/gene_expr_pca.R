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
  map(8) %>% # Change this later to match snakemake paths
  unlist()
colnames(means) <- expNames
print(means)

# Split into train/test
trainMeans <- means[,which(colnames(means) %in% trainTestSplit$train)]
testMeans <- means[,which(colnames(means) %in% trainTestSplit$test)]

trainMeans <- t(trainMeans)
testMeans <- t(testMeans)

ppcaTrain <- ppca(as.matrix(trainMeans), nPcs=20)
ppcaTest <- predict(ppcaTrain, newdata = as.matrix(testMeans))

trainExprScores <- ppcaTrain@scores
testExprScores <- ppcaTest$scores

saveRDS(list(trainExpr=trainExprScores, testExpr=testExprScores), "/home/campbell/cfang/automl_scrna/data/avg_expr_pca.RDS")