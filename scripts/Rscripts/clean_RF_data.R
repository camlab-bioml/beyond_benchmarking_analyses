suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(stringr)
  library(stringi)
})

args <- commandArgs(trailingOnly=TRUE)

designMat <- readRDS(args[[1]])
ari <- read.csv(args[[2]])
nclusts <- read.csv(args[[3]])
exprScores <- readRDS(args[[4]])

# Remove duplicated columns
designMat <- designMat[,-which(colnames(designMat) %in% c("sum.1", "detected.1", "total.1", "pct_counts_in_top_20_features", "percent.top_50", "total_counts","total_features"))]
designMat$pipelines <- gsub(pattern="\"", replacement="", designMat$pipelines)
#designMat <- merge(designMat, avg_expr_scores, by="name")
designMat <- designMat %>%
  dplyr::rename(filt = V2)%>%
  dplyr::rename(dims = V3)%>%
  dplyr::rename(norm = V4)%>%
  dplyr::rename(res = V5)

# Clean up ARI colnames
ari$X <- NULL
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.supervised\\.csv", replacement=""))
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.", replacement="-"))

# Merge number of clusters with design matrix
colnames(nclusts) <- unlist(lapply(colnames(nclusts), gsub, pattern="*\\.", replacement="-"))
nclusts$X <- gsub(pattern="\"", replacement="", nclusts$X)
nclusts$X <- gsub(pattern="\\;", replacement="\\.", nclusts$X)
nclusts$X <- gsub(pattern="\\=", replacement="\\.", nclusts$X)
colnames(nclusts)[which(colnames(nclusts)=="X")] <- "pipelines"

nclusts <- as_tibble(nclusts)
nclusts <- pivot_longer(nclusts, cols=which(grepl("^E-*", colnames(nclusts))))

designClusts <- merge(designMat, nclusts, by=c("pipelines", "name"))
colnames(designClusts)[which(colnames(designClusts)=="value")] <- "nclusts"
designMat <- designClusts

# Clean up expression PCA scores and merge with design matrix
avg_expr_scores <- rbind(exprScores$trainExpr, exprScores$testExpr)
avg_expr_scores <- avg_expr_scores[,1:20]
avg_expr_scores <- cbind(name=rownames(avg_expr_scores), avg_expr_scores)
rownames(avg_expr_scores) <- NULL
avg_expr_scores <- as_tibble(avg_expr_scores)
designMat <- merge(designMat, avg_expr_scores, by="name")

print(head(designMat))

saveRDS(designMat, "/home/campbell/cfang/automl_scrna/data/uncorrected_designMatrix_scaled_cleaned.RDS")
saveRDS(ari, "/home/campbell/cfang/automl_scrna/results/supervised_metrics/ari_unscaled_cleaned.RDS")
saveRDS(nclusts, "/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/num_clusters_cleaned.RDS")
