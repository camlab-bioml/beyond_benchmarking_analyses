suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
})
#snakemake rule inputs all metric files of one experiment at a time
expsMetrics <- commandArgs(trailingOnly = TRUE)
sil <- c()
db <- c()
ch <- c()
pipelines <- c()
print(expsMetrics)
print(length(expsMetrics))
# Each file represents one pipeline run
for (i in 1:length(expsMetrics)){
  metrics <- read.csv(expsMetrics[[i]])
  print(i)
  print(expsMetrics[[i]])
  print(metrics)
  sil <- c(sil, metrics$silhouette)
  db <- c(db, metrics$db_index)
  ch <- c(ch, metrics$ch_index)
  pipelines <- c(pipelines, metrics$pipeline)
}
print(sil)
print(pipelines)
metrics_matrix <- cbind(pipelines,sil,db,ch)

exp <- strsplit(expsMetrics, split="/metrics-")[[1]][[1]]
exp <- strsplit(exp, split="metrics/")[[1]][[2]]
print(metrics_matrix)

metricsFilename <- paste0("/home/campbell/cfang/automl_scrna/results/metrics/metrics_matrices/", exp)
metricsFilename <- paste0(metricsFilename, "-metrics.csv")
print(metricsFilename)
write.csv(metrics_matrix, metricsFilename)


