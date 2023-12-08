library(SingleCellExperiment)
library(umap)
library(Seurat)
library(scuttle)
library(tibble)
library(ggplot2)
library(ggforce)

# model <- readRDS("RFPreds_NumericParams.RDS")
# ari <- readRDS("ari_unscaled_cleaned.RDS")
# 
# 
# sil <- model$sil
# testMat <- cbind(sil$testData, preds=sil$testPreds)
# 
# bestDataset <- "E-MTAB-9221"
# best <- testMat[which(testMat$name==bestDataset),] %>%
#   select(name, pipelines, preds)
# 
# bestARIind <- which.max(best$preds)
# bestPip <- best$pipelines[which.max(best$preds)]
# 
# pipeAri <- ari[bestDataset]
# ari <- cbind(pipelines=ari$pipelines, ari=pipeAri)
# 
# which(ari$pipelines == bestPip)

bestPredLabels <- read.csv("E-MTAB-9221-bestPred-labelled.csv")
sce <- readRDS("E-MTAB-9221-SCE.RDS")
sce <- logNormCounts(sce)

set.seed(03112023)
sce <- scater::runPCA(sce, exprs_values="logcounts")
sce <- scater::runUMAP(sce, dimred="PCA")

umapLayout <- reducedDim(sce, "UMAP")


predLabs <- bestPredLabels[,3]
trueLabs <- bestPredLabels[,4]

trueLabs <- as.data.frame(cbind(cell=bestPredLabels$X, true=trueLabs))
predLabs <- as.data.frame(cbind(cell=bestPredLabels$X, pred=predLabs))

predLabs$pred <- as.character(as.numeric(predLabs$pred)+1)

umapLayout <- umapLayout %>%
  as.data.frame() %>%
  rownames_to_column(var="cell") %>% 
  merge(as.data.frame(predLabs), by="cell")%>%
  dplyr::rename(label=pred)

umapLayoutTrue <- merge(umapLayout, trueLabs, by="cell", all=TRUE) %>%
  select(cell, UMAP1, UMAP2, true)%>%
  dplyr::rename(label=true)

pred_ari_annot <- data.frame(UMAP1 = max(umapLayout$UMAP1)-3.5,
                                UMAP2 = max(umapLayout$UMAP2), 
                                ARI=ari[2,2])
pPred <- ggplot(umapLayout, aes(x=UMAP1, y=UMAP2, colour=label))+
  geom_point()+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("Best predicted pipeline")+
  geom_text(data=pred_ari_annot, aes(x=UMAP1, y=UMAP2, 
                                     label= paste("ARI:", round(ARI, 3))),
            color="black",
            size=5)+
  cowplot::theme_cowplot()

pTrue <- ggplot(umapLayoutTrue, aes(x=UMAP1, y=UMAP2, colour=label))+
  geom_point()+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("Ground truth labels")+
  cowplot::theme_cowplot()


set.seed(1132023)
randInt <- sample(1:324, 1)

randLabels <- read.csv("E-MTAB-9221-rand-labelled.csv")
randLabs <- as.data.frame(cbind(cell=randLabels$X, rand=randLabels[,3]))
randLabs$rand <- as.character(as.numeric(randLabs$rand)+1)

umapLayoutRand <-  merge(umapLayout, randLabs, by="cell", all=TRUE) %>%
  select(cell, UMAP1, UMAP2, rand)%>%
  dplyr::rename(label=rand)

rand_ari_annot <- data.frame(UMAP1 = max(umapLayout$UMAP1)-3.5,
                             UMAP2 = max(umapLayout$UMAP2), 
                             ARI=ari[312,2])
pRand <- ggplot(umapLayoutRand, aes(x=UMAP1, y=UMAP2, colour=label))+
  geom_point()+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("Randomly selected pipeline")+
  geom_text(data=rand_ari_annot, aes(x=UMAP1, y=UMAP2, 
                                     label= paste("ARI:", round(ARI, 3))),
            color="black",
            size=5)+
  cowplot::theme_cowplot()
  #facet_wrap(~ labelType)

pdf("bestPredVsRandomFewerGroundTruths.pdf", height=4, width=14.5)
pTrue + pPred + pRand
dev.off()

