suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
})
#snakemake rule inputs all metric files of one experiment at a time
args <- commandArgs(trailingOnly = TRUE)
metricsFilename <- args[[1]]
expsMetrics <- args[-1]

aris <- c()
mis <- c()
fms <- c()
vms <- c()
homos <- c()
comps <- c()
pipelines <- c()

print(expsMetrics)
print(length(expsMetrics))
# Each file represents one pipeline run
for (i in 1:length(expsMetrics)){
  metrics <- read.csv(expsMetrics[[i]])
  print(i)
  print(expsMetrics[[i]])
  print(metrics)
  aris <- c(aris, metrics$ARI)
  mis <- c(mis, metrics$mutual_info)
  homos <- c(homos, metrics$homogeneity)
  comps <- c(comps, metrics$completeness)
  vms <- c(vms, metrics$vmeasure)
  fms <- c(fms, metrics$FM)
  pipelines <- c(pipelines, metrics$run)
}
print(aris)
print(pipelines)
metrics_matrix <- cbind(pipelines,aris, mis, homos, comps, vms, fms)
print(metrics_matrix)

print(metricsFilename)
write.csv(metrics_matrix, metricsFilename)


