library(tidyverse)
library(here)

metadata <- read.csv(here("cleaned_collated_metadata.csv"))
metadata$X <- NULL

train <- readRDS(here("data", "corrected_RF_trainMatNumeric_unscaled.RDS"))
trainData <- train$data

trainWithMeta <- merge(trainData, metadata, by="name")

train$data <- trainWithMeta
saveRDS(train, here("data", "corrected_RF_trainMatNumeric_unscaled_withMeta.RDS"))

# test
test <- readRDS(here("data", "corrected_RF_testMatNumeric_unscaled.RDS"))
testData <- test$data

testWithMeta <- merge(testData, metadata, by="name")

test$data <- testWithMeta
saveRDS(test, here("data", "corrected_RF_testMatNumeric_unscaled_withMeta.RDS"))

design <- readRDS(here("data", "corrected_designMat_unscaled.RDS"))
designWithMeta <- merge(design, metadata, by="name")
saveRDS(designWithMeta, here("data", "corrected_designMat_unscaled_withMeta.RDS"))
