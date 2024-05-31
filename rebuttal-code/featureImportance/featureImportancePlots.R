library(randomForest)

silRF <- readRDS("data/original-models/RFTunedsil_NumericParams.RDS")
varImpPlot(silRF$finalModel, main="Variable Importance for SIL RF")


dbRF <- readRDS("data/original-models/RFTuneddb_NumericParams.RDS")
varImpPlot(dbRF$finalModel, main="Variable Importance for DB RF")

chRF <- readRDS("data/original-models/RFTunedch_NumericParams.RDS")
varImpPlot(chRF$finalModel, main="Variable Importance for CH RF")

gseaRF <- readRDS("data/original-models/RFTunedgsea_NumericParams.RDS")
varImpPlot(gseaRF$finalModel, main="Variable Importance for GSEA RF")

pdf("results/figures/randomForest/original-models/RFFeatImpPlots.pdf", height=15, width=10)
par(mfrow=c(2,2))
varImpPlot(chRF$finalModel, main="Variable Importance for CH RF")
varImpPlot(dbRF$finalModel, main="Variable Importance for DB RF")
varImpPlot(silRF$finalModel, main="Variable Importance for SIL RF")
varImpPlot(gseaRF$finalModel, main="Variable Importance for GSEA RF")

dev.off()
