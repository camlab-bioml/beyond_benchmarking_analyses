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

dataset <- "E-MTAB-8410"
best <- testMat[which(testMat$name==dataset),] %>%
  dplyr::select(name, pipelines, preds)

bestARIind <- which.max(best$preds)
bestPip <- best$pipelines[which.max(best$preds)]


pipeAri <- ari[dataset]
ari <- cbind(pipelines=ari$pipelines, ari=pipeAri)

which(ari$pipelines == bestPip)

bestPipPredictions <- read.csv("data/E-MTAB-8410-RF-best-predicted-pipeline-by-sil.csv")
seu <- readRDS("data/experiments/E-MTAB-8410-SEU.RDS")
sce <- as.SingleCellExperiment(seu)

#scater::plotUMAP(sce, colour_by="pct_counts_Mt") + scater::plotUMAP(sce, colour_by="seurat_clusters")


sce_pip <- sce[,bestPipPredictions$X]
stopifnot(all.equal(colnames(sce_pip), bestPipPredictions$X))
sce_pip$best_pip_preds <- as.factor(bestPipPredictions[,2])

sce$log10_detected <- log(sce$detected)
sce_pip$log10_detected <- log(sce_pip$detected)

author_labels <- read.table("data/ExpDesign-labels-E-MTAB-8410.tsv", sep="\t")
authlabs <- author_labels[,c("V1", "V22")]
colnames(authlabs) <- c("X", "author_labels")
authlabs <- authlabs[-1,]
authlabs <- authlabs[-which(authlabs$author_labels==""),]
authlabs$author_labels_num <- as.numeric(as.factor(authlabs$author_labels))

authlabs <- authlabs[authlabs$X %in% colnames(sce),]
sce_auth <- sce[,authlabs$X]
sce_auth$author_labels <- as.factor(authlabs$author_labels_num)
sce_auth$author_labels_name <- authlabs$author_labels



sce_auth_panel1 <- sce_auth[,sce_auth$author_labels==8]
p1_cell_inds <- colnames(sce_auth_panel1)
p1_cell_inds <- intersect(p1_cell_inds, colnames(sce_pip))

sce_pip_panel1 <- sce_pip[,p1_cell_inds]
sce_seu_panel1 <- sce[,p1_cell_inds]


seu_props <- table(sce_pip_panel1$seurat_clusters)/length(sce_pip_panel1$seurat_clusters)
pip_props <- table(sce_pip_panel1$best_pip_preds)/length(sce_pip_panel1$best_pip_preds)

seu_props <- seu_props[which(seu_props!=0)]
pip_props <- pip_props[which(pip_props!=0)]

seu_props <- as.data.frame(seu_props)
seu_props$method <- "Seurat"

pip_props <- as.data.frame(pip_props)
pip_props$method <- "Best predicted pipeline"

all_props <- rbind(seu_props, pip_props)

barp2 <- ggplot(all_props, aes(y=Freq, x=method, fill=Var1))+
  geom_bar(position="fill", stat="identity")+
  ylab("Cluster Proportions")+
  xlab("Method")+
  cowplot::theme_cowplot()+
  theme(legend.position="none")+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        title=element_text(hjust=0.5, size=25))+
  ggtitle("Consensus molecular\nsubtype 3")




#df <- data.frame(seu_props, pip_props)

## panel 2 ###
sce_auth_panel2 <- sce_auth[,sce_auth$author_labels %in% c(40, 15)]
p2_cell_inds <- colnames(sce_auth_panel2)
p2_cell_inds <- intersect(p2_cell_inds, colnames(sce_pip))

sce_pip_panel2 <- sce_pip[,p2_cell_inds]
sce_seu_panel2 <- sce[,p2_cell_inds]



seu_props <- table(sce_pip_panel2$seurat_clusters)/length(sce_pip_panel2$seurat_clusters)
pip_props <- table(sce_pip_panel2$best_pip_preds)/length(sce_pip_panel2$best_pip_preds)

seu_props <- seu_props[which(seu_props!=0)]
pip_props <- pip_props[which(pip_props!=0)]

seu_props <- as.data.frame(seu_props)
seu_props$method <- "Seurat"

pip_props <- as.data.frame(pip_props)
pip_props$method <- "Best predicted pipeline"

all_props <- rbind(seu_props, pip_props)

barp3 <- ggplot(all_props, aes(y=Freq, x=method, fill=Var1))+
  geom_bar(position="fill", stat="identity")+
  ylab("Cluster Proportions")+
  xlab("Method")+
  cowplot::theme_cowplot()+
  theme(legend.position="none")+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        title=element_text(hjust=0.5, size=25))+
  ggtitle("IgA+ plasma cells")


#### panel 3 ####
sce_auth_panel3 <- sce_auth[,sce_auth$author_labels %in% c(32)]
p3_cell_inds <- colnames(sce_auth_panel3)
p3_cell_inds <- intersect(p3_cell_inds, colnames(sce_pip))

sce_pip_panel3 <- sce_pip[,p3_cell_inds]
sce_seu_panel3 <- sce[,p3_cell_inds]


seu_props <- table(sce_pip_panel3$seurat_clusters)/length(sce_pip_panel3$seurat_clusters)
pip_props <- table(sce_pip_panel3$best_pip_preds)/length(sce_pip_panel3$best_pip_preds)

seu_props <- seu_props[which(seu_props!=0)]
pip_props <- pip_props[which(pip_props!=0)]

seu_props <- as.data.frame(seu_props)
seu_props$method <- "Seurat"

pip_props <- as.data.frame(pip_props)
pip_props$method <- "Best predicted pipeline"

all_props <- rbind(seu_props, pip_props)

barp4 <- ggplot(all_props, aes(y=Freq, x=method, fill=Var1))+
  geom_bar(position="fill", stat="identity")+
  ylab("Cluster Proportions")+
  xlab("Method")+
  cowplot::theme_cowplot()+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        title=element_text(hjust=0.5, size=25))+
  ggtitle("Stromal 3")

pdf("results/figures/E-MTAB-8410-default-seurat-comparison-barplots.pdf", height=8, width=15)
barp2 + barp3 + barp4
dev.off()

# ggsave("results/figures/E-MTAB-8410-best_pip_sil_vs_seurat_default_panel1.png",
#        height=4, width=13.5, dpi=400)
#
# scater::plotUMAP(sce_pip, colour_by="best_pip_preds", text_by="best_pip_preds")+ scater::plotUMAP(sce, colour_by="seurat_clusters", text_by="seurat_clusters") +scater::plotUMAP(sce_auth, text_by="author_labels_name", colour_by="author_labels")+labs(colour="Author labels")# scale_colour_manual(values=pal)
# dev.off()
#
# ggsave("results/figures/E-MTAB-8410-best_pip_sil_vs_seurat_default_panel2.png",
#        height=4, width=13.5, dpi=400)
# scater::plotUMAP(sce_pip_panel1, colour_by="best_pip_preds", text_by="best_pip_preds")+labs(colour="Best predicted pipeline labels")+ scater::plotUMAP(sce_seu_panel1, colour_by="seurat_clusters", text_by="seurat_clusters")+labs(colour="Default Seurat labels") +scater::plotUMAP(sce_auth_panel1, colour_by="author_labels_name", text_by="author_labels_name")+labs(colour="Author labels")# scale_colour_manual(values=pal)
# dev.off()
#
#
# ggsave("results/figures/E-MTAB-8410-best_pip_sil_vs_seurat_default_panel3.png",
#        height=4, width=13.5, dpi=400)
# scater::plotUMAP(sce_pip_panel2, colour_by="best_pip_preds", text_by="best_pip_preds")+ scater::plotUMAP(sce_seu_panel2, colour_by="seurat_clusters", text_by="seurat_clusters") +scater::plotUMAP(sce_auth_panel2, colour_by="author_labels_name", text_by="author_labels_name")+labs(colour="Author labels")# scale_colour_manual(values=pal)
# dev.off()
#
# ggsave("results/figures/E-MTAB-8410-best_pip_sil_vs_seurat_default_panel4.png",
#        height=4, width=13.5, dpi=400)
# scater::plotUMAP(sce_pip_panel3, colour_by="best_pip_preds", text_by="best_pip_preds")+labs(colour="Best predicted pipeline labels")+ scater::plotUMAP(sce_seu_panel3, colour_by="seurat_clusters", text_by="seurat_clusters")+labs(colour="Default Seurat labels") +scater::plotUMAP(sce_auth_panel3, colour_by="author_labels_name", text_by="author_labels_name")+labs(colour="Author labels")
# dev.off()
#
#
#
