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

# Get average expression for each dataset
means <- args %>%
  map(readRDS) %>%
  map(counts) %>%
  map(rowMeans) %>%
  map(as.data.frame)%>%
  map(rownames_to_column,var="gene") %>%
  purrr::reduce(merge, by="gene", all=TRUE)%>%
  column_to_rownames(var="gene")

# Get experiment names
expNames <- strsplit(args, split="/") %>%
  map(3) %>% # Change this later to match snakemake paths
  unlist()
colnames(means) <- expNames
#print(means)

#print(means)

means <- t(means)
means[is.na(means)] <- 0

print(head(means))

print(head(colnames(means)))
pcaTrain <- readRDS("/home/campbell/cfang/bb-rebuttal/results/gene_expr_pca_model.RDS")
print(pcaTrain)
loadings <- pcaTrain@loadings
loadings <- loadings[rownames(loadings) %in% colnames(means),]

means <- means[,colnames(means) %in% rownames(loadings)]
print('means')
print(dim(as.matrix(means)))
print('loadings')
print(dim(loadings))

testExprScores <- as.matrix(means) %*% loadings

print(testExprScores)


saveRDS(list(testExpr=testExprScores), "/home/campbell/cfang/bb-rebuttal/data/large_samples_avg_expr_pca.RDS")
