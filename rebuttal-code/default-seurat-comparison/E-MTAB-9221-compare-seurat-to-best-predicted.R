library(SingleCellExperiment)
library(umap)
library(Seurat)
library(scuttle)
library(tibble)
library(ggplot2)
library(ggforce)
library(here)
library(dplyr)

model <- readRDS(here("data", "RFPreds_NumericParams.RDS"))
ari <- readRDS(here("data", "ari_unscaled_cleaned.RDS"))


sil <- model$sil
testMat <- cbind(sil$testData, preds=sil$testPreds)

dataset <- "E-MTAB-9221"
best <- testMat[which(testMat$name==dataset),] %>%
  dplyr::select(name, pipelines, preds)

bestARIind <- which.max(best$preds)
bestPip <- best$pipelines[which.max(best$preds)]


pipeAri <- ari[dataset]
ari <- cbind(pipelines=ari$pipelines, ari=pipeAri)

which(ari$pipelines == bestPip)

bestPipPredictions <- read.csv("data/E-MTAB-9221-RF-best-predicted-pipeline-by-sil.csv")
seu <- readRDS("data/experiments/E-MTAB-9221-SEU.RDS")
sce <- as.SingleCellExperiment(seu)

#scater::plotUMAP(sce, colour_by="pct_counts_Mt") + scater::plotUMAP(sce, colour_by="seurat_clusters")


sce_pip <- sce[,bestPipPredictions$X]
stopifnot(all.equal(colnames(sce_pip), bestPipPredictions$X))
sce_pip$best_pip_preds <- as.factor(bestPipPredictions[,2])

sce$log10_detected <- log(sce$detected)
sce_pip$log10_detected <- log(sce_pip$detected)

author_labels <- read.table("data/ExpDesign-labels-E-MTAB-9221.tsv", sep="\t")
authlabs <- author_labels[,c("V1", "V24")]
colnames(authlabs) <- c("X", "author_labels")
authlabs <- authlabs[-1,]
authlabs <- authlabs[-which(authlabs$author_labels==""),]
authlabs$author_labels_num <- as.numeric(as.factor(authlabs$author_labels))

authlabs <- authlabs[authlabs$X %in% colnames(sce),]
sce_auth <- sce[,authlabs$X]
sce_auth$author_labels <- as.factor(authlabs$author_labels_num)
sce_auth$author_labels_name <- authlabs$author_labels


ggsave("results/figures/E-MTAB-9221-best_pip_sil_vs_seurat_default_panel1.png",
       height=4, width=13.5, dpi=400)

scater::plotUMAP(sce_pip, colour_by="best_pip_preds", text_by="best_pip_preds")+ scater::plotUMAP(sce, colour_by="seurat_clusters", text_by="seurat_clusters") +scater::plotUMAP(sce_auth, text_by="author_labels_name", colour_by="author_labels_name")+labs(colour="Author labels")# scale_colour_manual(values=pal)
dev.off()






