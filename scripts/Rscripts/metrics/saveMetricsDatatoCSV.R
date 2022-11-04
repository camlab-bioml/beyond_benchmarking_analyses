suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(DT) # CRAN
  library(Matrix) # CRAN
  library(Metrics)
  library(pipeComp)
})

args <- commandArgs(trailingOnly=TRUE)

numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}

numClusters <- read.csv(args[[1]]) # read in results/pipecomp_outputs/num_clusters.csv
pipelines <- parsePipNames(numClusters$X)
pipelines <- pipelines[-which(pipelines$resolution == 0.01),]
rownames(pipelines) <- numClusters$X[-grep("0.01", numClusters$X)]
rownames(pipelines) <- gsub(pattern="\"", replacement="", rownames(pipelines))
rownames(pipelines) <- gsub(pattern="\\;", replacement="\\.", rownames(pipelines))
rownames(pipelines) <- gsub(pattern="\\=", replacement="\\.", rownames(pipelines))

print(head(pipelines))
write.csv(pipelines, file="results/paperData/pipelineParams.csv")

colData <- readRDS(args[[2]]) # read in data/experiments/colData.RDS

colData %>% 
  as_tibble() %>%
  mutate_if(numbers_only, as.numeric)

colnames(colData)[[1]] <- "name"
print(head(colData))

exprScores <- readRDS(args[[8]])
avg_expr_scores <- rbind(exprScores$trainExpr, exprScores$testExpr)
avg_expr_scores <- avg_expr_scores[,1:20]
avg_expr_scores <- cbind(name=rownames(avg_expr_scores), avg_expr_scores)
rownames(avg_expr_scores) <- NULL
avg_expr_scores <- as_tibble(avg_expr_scores)

colData <- colData[,-which(duplicated(t(colData)))]
colData <- merge(colData, avg_expr_scores, by="name")

for (i in 1:length(colnames(colData))){
  if (grepl("^V", colnames(colData)[[i]])){
    print(colnames(colData))[[i]]
    num <- as.numeric(strsplit(colnames(colData)[[i]], split="V")[[1]][[2]])-1
    colnames(colData)[[i]] <- paste0("PC", num)
  }
}
print(head(colData))
write.csv(colData, file="results/paperData/datasetFeatures.csv") #add gene expr pca?

silUnscaled <- read.csv(args[[3]]) # read in results/metrics/sil_unscaled.csv
silUnscaled <- silUnscaled[-grep("0.01", silUnscaled$pipelines),]

dbUnscaled <- read.csv(args[[4]]) # read in results/metrics/db_unscaled.csv
dbUnscaled <- dbUnscaled[-grep("0.01", dbUnscaled$pipelines),]

chUnscaled <- read.csv(args[[5]]) # read in results/metrics/ch_unscaled.csv
print(head(chUnscaled))
chUnscaled <- chUnscaled[-grep("0.01", chUnscaled$pipelines),]

gseaUnscaled <- read.csv(args[[6]]) 
gseaUnscaled$X <- NULL

gseaUnscaled <- gseaUnscaled[-grep("0.01", gseaUnscaled$pipelines),]
gseaUnscaled$pipelines <- gsub(pattern="\"", replacement="", gseaUnscaled$pipelines)
gseaUnscaled$pipelines <- gsub(pattern="\\;", replacement="\\.", gseaUnscaled$pipelines)
gseaUnscaled$pipelines <- gsub(pattern="\\=", replacement="\\.", gseaUnscaled$pipelines)
print(head(gseaUnscaled))
gseaUnscaled <- gseaUnscaled %>%
  pivot_wider(names_from=name, values_from=means)

print(head(gseaUnscaled))
write.csv(silUnscaled, "results/paperData/sil_unscaled.csv")
write.csv(dbUnscaled, "results/paperData/db_unscaled.csv")
write.csv(chUnscaled, "results/paperData/ch_unscaled.csv")
write.csv(gseaUnscaled, "results/paperData/gsea_unscaled.csv")

designMat <- readRDS(args[[7]]) # read in "corrected_designMat_unscaled.RDS"

sil <- designMat %>%
  select(pipelines, name, silCorImp) %>%
  pivot_wider(names_from=name, values_from=silCorImp)

db <- designMat %>%
  select(pipelines, name, dbCorImp) %>%
  pivot_wider(names_from=name, values_from=dbCorImp)

ch <- designMat %>%
  select(pipelines, name, chCorImp) %>%
  pivot_wider(names_from=name, values_from=chCorImp)

gsea <- designMat %>%
  select(pipelines, name, gseaScaledImp) %>%
  pivot_wider(names_from=name, values_from=gseaScaledImp)

write.csv(sil, "results/paperData/silCorrectedImputed.csv")
write.csv(db, "results/paperData/dbCorrectedImputed.csv")
write.csv(ch, "results/paperData/chCorrectedImputed.csv")
write.csv(gsea, "results/paperData/gseaScaledImputed.csv")

