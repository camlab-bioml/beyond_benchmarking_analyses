suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(ggplot2)
  library(stringr)
  library(stringi)
  library(ggforce)
})

args <- commandArgs(trailingOnly=TRUE)
print(args)

#RFmodels <- readRDS(args[[1]])
LRmodels <- readRDS(args[[5]])

silRF <- readRDS(args[[1]])
dbRF <- readRDS(args[[2]])
chRF <- readRDS(args[[3]])
gseaRF <- readRDS(args[[4]])


print("RF")
print(names(silRF))
pdf("results/figures/randomForest/RFTuning.pdf")
plot(silRF)
plot(dbRF)
plot(chRF)
plot(gseaRF)
dev.off()

silLR <- LRmodels$sil$model
dbLR <- LRmodels$db$model
chLR <- LRmodels$ch$model
gseaLR <- LRmodels$gsea$model

print("LR")
pdf("results/figures/linearRegression/LRTuning.pdf")
plot(silLR)
plot(dbLR)
plot(chLR)
plot(gseaLR)
dev.off()