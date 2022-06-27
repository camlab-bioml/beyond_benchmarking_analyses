library(caret)
library(tidyverse)
#models <- readRDS("RFPreds_NumericParams.RDS")
args <- commandArgs(trailingOnly = TRUE)
models <- readRDS(args[[1]])

sil <- varImp(models$sil$model, scale=TRUE)
db <- varImp(models$db$model, scale=TRUE)
ch <- varImp(models$ch$model, scale=TRUE)
gsea <- varImp(models$gsea$model, scale=TRUE)

silRanks <- rownames(sil)[order(-sil)] %>%
  as_tibble()%>%
  mutate(rank=row_number())
dbRanks <- rownames(db)[order(-db)]%>%
  as_tibble()%>%
  mutate(rank=row_number())
chRanks <- rownames(ch)[order(-ch)]%>%
  as_tibble()%>%
  mutate(rank=row_number())
gseaRanks <- rownames(gsea)[order(-gsea)]%>%
  as_tibble()%>%
  mutate(rank=row_number())


allImp <- cbind(sil, db, ch, gsea)
colnames(allImp) <- c("sil", "db", "ch", "gsea")
topTen <- rowMeans(allImp)[order(-rowMeans(allImp))]
topTen <- names(topTen)[1:10]

silTopTen <- silRanks %>%
  filter(value %in% topTen)
dbTopTen <- dbRanks %>%
  filter(value %in% topTen)
chTopTen <- chRanks %>%
  filter(value %in% topTen)
gseaTopTen <- gseaRanks %>%
  filter(value %in% topTen)

ranksDF <- merge(silTopTen, dbTopTen, by="value")
ranksDF <- merge(ranksDF, chTopTen, by="value")
ranksDF <- merge(ranksDF, gseaTopTen, by="value")
ranksDF <- column_to_rownames(ranksDF, "value")
colnames(ranksDF) <- c("SIL", "DB", "CH", "GSEA")


for (i in 1:length(rownames(ranksDF))){
  if (grepl("^V", rownames(ranksDF)[[i]])){
    print(rownames(ranksDF))[[i]]
    num <- as.numeric(strsplit(rownames(ranksDF)[[i]], split="V")[[1]][[2]])-1
    rownames(ranksDF)[[i]] <- paste0("PC", num)
  }
}

print(ranksDF)
f1 = circlize::colorRamp2(seq(1, 10, length = nrow(ranksDF)), hcl.colors(nrow(ranksDF),"viridis"))

pdf("/home/campbell/cfang/automl_scrna/results/figures/RFPlots/RF_featImp.pdf")
#pdf("PaperFigures/RF_featImp.pdf")
ComplexHeatmap::Heatmap(as.matrix(ranksDF), col=f1,row_dend_reorder = TRUE, 
                        name="Importance\nrank",
                        row_title="Feature",
                        column_title="Metric",
                        column_title_side = "bottom",
                        heatmap_legend_param = list(at = c(2,4,6,8,10)
                        ))
dev.off()