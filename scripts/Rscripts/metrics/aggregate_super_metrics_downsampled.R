suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(scran) # BioConductor
  library(purrr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
downsampleType <- args[[1]]
expsMatrices <- args[-1]

ari <- c()
mi <- c()
homo <- c()
fm <- c()
vm <- c()
comp <- c()

for (i in 1:length(expsMatrices)){
  print(expsMatrices[[i]])
  metrics <- read.csv(expsMatrices[[i]])
  colnames(metrics)[[1]] <- "pipelines"
  expName <- strsplit(expsMatrices[[i]], "-metrics.csv")[[1]][[1]]
  expName <- strsplit(expName, "metrics_matrices/")[[1]][[2]]
  #print(expName)
  #print(dim(metrics))
  names(metrics)[[1]] <- "index"
  #print(head(metrics))
  print(any(str_length(metrics$pipelines) < 10))
  print(metrics[which(str_length(metrics$pipelines) < 10),])
  
  expNames <- rep(expName, nrow(metrics))
  
  metrics$pipelines <- gsub(pattern="\\+", replacement="\\.", metrics$pipelines)
  metrics$pipelines <- gsub(pattern="\\=", replacement = "\\.", metrics$pipelines)
  
  ari[[i]] <- data.frame("pipelines"=metrics$pipelines, "ari"=metrics$aris, "exp"=expNames)
  mi[[i]] <- data.frame("pipelines"=metrics$pipelines, "mi"=metrics$mis, "exp"=expNames)
  homo[[i]] <- data.frame("pipelines"=metrics$pipelines, "homo"=metrics$homos, "exp"=expNames)
  fm[[i]] <- data.frame("pipelines"=metrics$pipelines, "fm"=metrics$fms, "exp"=expNames)
  vm[[i]] <- data.frame("pipelines"=metrics$pipelines, "vm"=metrics$vms, "exp"=expNames)
  comp[[i]] <- data.frame("pipelines"=metrics$pipelines, "comp"=metrics$comp, "exp"=expNames)
}
# do.call(rbind(sil)), then make it a tibble and pivot
print(unique(ari$pipelines))
ariMat <- as_tibble(do.call(rbind, ari))
miMat <- as_tibble(do.call(rbind, mi))
homoMat <- as_tibble(do.call(rbind, homo))
fmMat <- as_tibble(do.call(rbind, fm))
vmMat <- as_tibble(do.call(rbind, vm))
compMat <- as_tibble(do.call(rbind, comp))

ariMat <- pivot_wider(ariMat, names_from = "exp", values_from="ari")
miMat <- pivot_wider(miMat, names_from = "exp", values_from="mi")
homoMat <- pivot_wider(homoMat, names_from = "exp", values_from="homo")
fmMat <- pivot_wider(fmMat, names_from = "exp", values_from="fm")
vmMat <- pivot_wider(vmMat, names_from = "exp", values_from="vm")
compMat <- pivot_wider(compMat, names_from = "exp", values_from="comp")

ariName <- "/home/campbell/cfang/automl_scrna/results/supervised_metrics/ari_unscaled_with_downsampled.csv"
miName <- "/home/campbell/cfang/automl_scrna/results/supervised_metrics/mi_unscaled_with_downsampled.csv"
homoName <- "/home/campbell/cfang/automl_scrna/results/supervised_metrics/homo_unscaled_with_downsampled.csv"
fmName <- "/home/campbell/cfang/automl_scrna/results/supervised_metrics/fm_unscaled_with_downsampled.csv"
vmName <- "/home/campbell/cfang/automl_scrna/results/supervised_metrics/vm_unscaled_with_downsampled.csv"
compName <- "/home/campbell/cfang/automl_scrna/results/supervised_metrics/comp_unscaled_with_downsampled.csv"

write.csv(ariMat, ariName)
write.csv(miMat, miName)
write.csv(homoMat, homoName)
write.csv(fmMat, fmName)
write.csv(vmMat, vmName)
write.csv(compMat, compName)


