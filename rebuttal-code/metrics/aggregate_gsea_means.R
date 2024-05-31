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
gseas <- args
data <- gseas %>%
  map(read_csv) %>%
  bind_rows() %>%
  select(pipelines, name, means)
print(data)
write.csv(data, "/home/campbell/cfang/bb-rebuttal/results/gsea_results/large_samples_all_gsea_means.csv")
