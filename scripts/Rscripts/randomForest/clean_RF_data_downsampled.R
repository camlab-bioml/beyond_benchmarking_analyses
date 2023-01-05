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
nclustsDownsample <- read.csv(args[[4]])
exprScores <- readRDS(args[[5]])

# Remove duplicated columns
print(colnames(designMat))
designMat <- designMat[,-which(colnames(designMat) %in% c("sum.1", "detected.1", "total.1", "pct_counts_in_top_20_features", "percent.top_50", "total_counts","total_features"))]
designMat$pipelines <- gsub(pattern="\"", replacement="", designMat$pipelines)
#designMat <- merge(designMat, avg_expr_scores, by="name")
designMat <- designMat %>%
  dplyr::rename(filt = V2)%>%
  dplyr::rename(dims = V3)%>%
  dplyr::rename(norm = V4)%>%
  dplyr::rename(res = V5)

#head(designMat)

# concatenate experiment ID and downsampleType
designMat <- designMat %>%
  mutate(name = paste(name, downsampleType, sep="-"))

#print(head(designMat))

# Clean up ARI colnames
ari$X <- NULL
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.supervised\\.csv", replacement=""))
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.", replacement="-"))

# Merge number of clusters for regular and downsampled experiments
nclusts <- merge(nclusts, nclustsDownsample, by="X")
  
# Merge number of clusters with design matrix
colnames(nclusts) <- unlist(lapply(colnames(nclusts), gsub, pattern="*\\.", replacement="-"))
nclusts$X <- gsub(pattern="\"", replacement="", nclusts$X)
nclusts$X <- gsub(pattern="\\;", replacement="\\.", nclusts$X)
nclusts$X <- gsub(pattern="\\=", replacement="\\.", nclusts$X)
colnames(nclusts)[which(colnames(nclusts)=="X")] <- "pipelines"

#head(nclusts)

nclusts <- as_tibble(nclusts)
colnames(nclusts) <-gsub(pattern="\\.", replacement="-", colnames(nclusts))
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
avg_expr_scores$name <- gsub(pattern = "[.]-SCE$" ,replacement = "", x  = avg_expr_scores$name)

designMat <- merge(designMat, avg_expr_scores, by="name")

# Remove duplicated predictors
designMat$total <- NULL
designMat$pct_coding <- NULL
designMat$pct_ribosomal <- NULL
designMat$pct_counts_Mt <- NULL

print(head(designMat))
# 
# saveRDS(designMat, "/home/campbell/cfang/automl_scrna/data/uncorrected_designMatrix_scaled_cleaned_with_downsampled.RDS")
# saveRDS(ari, "/home/campbell/cfang/automl_scrna/results/supervised_metrics/ari_unscaled_cleaned_with_downsampled.RDS")
# saveRDS(nclusts, "/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/num_clusters_cleaned.RDS_with_downsampled")
