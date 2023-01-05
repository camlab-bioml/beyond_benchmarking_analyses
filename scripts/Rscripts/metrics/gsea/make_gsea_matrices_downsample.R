suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(stringr)
  library(stringi)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
output <-args[[1]]
gseas <- args[-c(1,length(args))]

print(args)
print(length(args))
data <- gseas %>%
  map(read_csv) %>%
  bind_rows()
print(data)

print(tail(args, 1))
clusters <- readRDS(tail(args,1))
data$pipelines <- names(clusters)


exp.name <- strsplit(output, split="/")[[1]][[4]]
exp.name <- strsplit(exp.name, split="-gsea.csv")[[1]][[1]]
data$name <- rep(exp.name, 324)
print(exp.name)
#file.name <- paste0("/home/campbell/cfang/automl_scrna/results/gsea_results/gsea_matrices/", exp.name)
#file.name <- paste0(file.name, "-gsea.csv")
#print(file.name)
write.csv(data, output)