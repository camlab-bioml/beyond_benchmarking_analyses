suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
})
#snakemake rule inputs all metric files of one experiment at a time
args <- commandArgs(trailingOnly = TRUE)
output <- args[[1]]
expsMetrics <- args[-1]
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

print(metrics_matrix)

write.csv(metrics_matrix, output)


