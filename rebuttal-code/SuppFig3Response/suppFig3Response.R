library(tidyverse)
library(here)


# In supplementary Figure 3 - despite the lack of consensus over data sets,
# are the highest-ranking pipelines on one metric consistent between metrics?


metrics <- readRDS(here("data", "corrected_designMat_unscaled.RDS")) %>%
  select(silCorImp, dbCorImp, chCorImp, gseaScaledImp, pipelines)

silBestPip <- metrics$pipelines[which.max(metrics$silCorImp)]
dbBestPip <- metrics$pipelines[which.max(metrics$dbCorImp)]
chBestPip <- metrics$pipelines[which.max(metrics$chCorImp)]
gseaBestPip <- metrics$pipelines[which.max(metrics$gseaScaledImp)]

silBestPip == dbBestPip
silBestPip == chBestPip
chBestPip == gseaBestPip
gseaBestPip == silBestPip
