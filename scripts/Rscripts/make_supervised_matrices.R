suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
})
#snakemake rule inputs all metric files of one experiment at a time
expsMetrics <- commandArgs(trailingOnly = TRUE)

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

exp <- strsplit(expsMetrics, split="metrics/")[[1]][[2]]
exp <- strsplit(exp, split="/")[[1]][[1]]
print(exp)
print(metrics_matrix)

metricsFilename <- paste0("/home/campbell/cfang/automl_scrna/results/supervised_metrics/metrics_matrices/", exp)
metricsFilename <- paste0(metricsFilename, "-supervised.csv")
print(metricsFilename)
write.csv(metrics_matrix, metricsFilename)


