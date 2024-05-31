library(tidyverse)
library(here)
metadata <- read.csv(here("collated_metadata_filledin.csv"), na.strings=c("", NA))


metadata_test <- metadata %>%
  select(name, grep("isolation", colnames(metadata)),
         grep("FromBlood", colnames(metadata)),
         grep("bias", colnames(metadata)))%>%
  mutate(isolation=coalesce(Comment.single.cell.isolation., Comment..single.cell.isolation.))%>%
  mutate(endBias=coalesce(Comment.end.bias., Comment..end.bias.))%>%
  select(c("name", "isolation" , "FromBlood", "endBias"))

metadata_test <- metadata_test %>%
  mutate(isolation=ifelse(grepl("10", metadata_test$isolation), "10X", isolation))%>%
  mutate(isolation=ifelse(grepl("FACS", metadata_test$isolation) | grepl("fluorescence", metadata_test$isolation), "FACS", isolation))%>%
  mutate(endBias=ifelse(grepl("3", metadata_test$endBias), "3prime", endBias))%>%
  mutate(endBias=ifelse(grepl("5", metadata_test$endBias), "5prime", endBias))%>%
  mutate(endBias=ifelse(grepl("none", metadata_test$endBias), "none", endBias))%>%
  mutate(isolation=ifelse(grepl("magnetic", metadata_test$isolation), "MACS", isolation))

# fill in some missing info
metadata_test$name[which(is.na(metadata_test$isolation))]
# [1] "E-GEOD-36552" "E-MTAB-6308"
metadata_test$name[which(is.na(metadata_test$endBias))]
# [1] "E-GEOD-36552" "E-MTAB-6308"

metadata_test$isolation[which(metadata_test$name=="E-GEOD-36552")] <- "manual"
metadata_test$endBias[which(metadata_test$name=="E-GEOD-36552")] <- "3prime"


metadata_test$isolation[which(metadata_test$name=="E-MTAB-6308")] <- "MACS"
metadata_test$endBias[which(metadata_test$name=="E-MTAB-6308")] <- "3prime"

metadata_test <- metadata_test %>%
  mutate(isFACS=ifelse(isolation=="FACS", "yes", "no"))%>%
  select(c(name, FromBlood, endBias, isFACS))


write.csv(metadata_test, "cleaned_collated_metadata.csv")
