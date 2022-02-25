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
write.csv(data, "/home/campbell/cfang/automl_scrna/results/gsea_results/all_gsea_means.csv")