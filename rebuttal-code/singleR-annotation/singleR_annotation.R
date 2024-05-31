library(SingleR)
library(here)
library(rhdf5)
library(zellkonverter)
#
# RFPreds_NumericParams <- readRDS("~/Documents/2022_autoMLscRNA_Cindy/data/RFPreds_NumericParams.RDS")
# testDatasets <- unique(RFPreds_NumericParams$sil$testData$name)
# rm(RFPreds_NumericParams)


sce <- readRDS(here("data", "experiments", "E-GEOD-114530-SCE.RDS"))
sce <- scuttle::logNormCounts(sce)

ref.sce <- readH5AD("data/Fetal_full_v3.h5ad")
ref.sce <- scuttle::logNormCounts(ref.sce)
rownames(ref.sce) <- rowData(ref.sce)$ID


preds<- SingleR(test=sce, ref=ref.sce, labels=ref.sce$celltype, de.method="wilcox")
table(preds$labels)
saveRDS(preds, here("data", "E-GEOD-114530-singleR-results.RDS"))
saveRDS(ref.sce, "data/HCA-fetal-kidney-atlas.RDS")
