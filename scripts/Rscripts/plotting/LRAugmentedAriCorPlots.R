library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)


preds <- readRDS("linearRegressionPreds_withDownsampled.RDS")
oldPreds <- readRDS("linearRegressionPreds.RDS")

ari <- read.csv("ari_unscaled_with_downsampled.csv")
oldAri <- readRDS("ari_unscaled_cleaned.RDS")


ari <- ari %>%
  rename_with(~ gsub(".supervised.csv", "", .x))

colnames(ari)[3:ncol(ari)] <- unlist(lapply(lapply(strsplit(colnames(ari)[3:ncol(ari)], split="\\."), "[", 4:9), FUN=paste, collapse="-"))
colnames(ari) <- gsub("-NA", "", colnames(ari))
ari$X <- NULL

allAri <- merge(ari, oldAri, by="pipelines")
allAri <- allAri %>%
  pivot_longer(cols=starts_with("E"), names_to="name")

sil <- preds$sil
db <- preds$db
ch <- preds$ch
gsea <- preds$gsea

silOrig <- oldPreds$sil
dbOrig <- oldPreds$db
chOrig <- oldPreds$ch
gseaOrig <- oldPreds$gsea

ariCor <- function(metric, ari){
  x <- cbind(metric$testData, preds=metric$testPreds)
  x <- merge(x, ari, by=c("pipelines", "name"))
  cor <- x %>%
    filter(!grepl("SCE", .$name)) %>%
    filter(!grepl("downsampleCounts", .$name)) %>%
    group_by(name) %>%
    group_modify(~ as.data.frame(cor(as.numeric(.$preds), as.numeric(.$value))))
  return(cor)
}

silCor <- ariCor(sil, allAri)
dbCor <- ariCor(db, allAri)
chCor <- ariCor(ch, allAri)
gseaCor <- ariCor(gsea, allAri)

colnames(silCor)[[2]] <- "sil"
colnames(dbCor)[[2]] <- "db"
colnames(chCor)[[2]] <- "ch"
colnames(gseaCor)[[2]] <- "gsea"

silCorOrig <- ariCor(silOrig, allAri)
dbCorOrig <- ariCor(dbOrig, allAri)
chCorOrig <- ariCor(chOrig, allAri)
gseaCorOrig <- ariCor(gseaOrig, allAri)

colnames(silCorOrig)[[2]] <- "sil"
colnames(dbCorOrig)[[2]] <- "db"
colnames(chCorOrig)[[2]] <- "ch"
colnames(gseaCorOrig)[[2]] <- "gsea"


allCor <- merge(silCor, dbCor, by="name")
allCor <- merge(allCor, chCor, by="name")
allCor <- merge(allCor, gseaCor, by="name")

allCorOrig <- merge(silCorOrig, dbCorOrig, by="name")
allCorOrig <- merge(allCorOrig, chCorOrig, by="name")
allCorOrig <- merge(allCorOrig, gseaCorOrig, by="name")

allCorOrig$modelType <- "original"
allCor$modelType <- "augmented"

allCor <- pivot_longer(allCor, cols =c(sil,db,ch,gsea) , names_to="metric")
allCorOrig <- pivot_longer(allCorOrig, cols =c(sil,db,ch,gsea) , names_to="metric")

bothModelsCor <- rbind(allCor, allCorOrig)

allCorSig <- allCor %>%
  group_by(metric, modelType)%>%
  group_modify(~as.data.frame(cbind(wilcox.test(x=as.numeric(.$value),  alternative="greater")$p.value)))

allCorSigOrig <- allCorOrig %>%
  group_by(metric, modelType)%>%
  group_modify(~as.data.frame(cbind(wilcox.test(x=as.numeric(.$value),  alternative="greater")$p.value)))

allCorPvals <- rbind(allCorSig, allCorSigOrig)


ggplot(bothModelsCor, aes(interaction(modelType, metric), value, fill=metric))+
  geom_boxplot()+
  ylab("correlation")+
  cowplot::theme_cowplot()+
  theme(axis.text.x=element_text(angle=45))+
  ggtitle("Correlation between LR predictions and ARI on test set")+
  geom_text(
    data=allCorPvals,
    mapping = aes(interaction(modelType, metric), y=0.75),
    label=round(allCorPvals$V1,10)
  )


# ggplot(merge(cbind(gsea$testData, preds=gsea$testPreds), allAri), aes(x=value, y=as.numeric(gseaScaledImp)))+
#   geom_point()+
#   facet_wrap(~name)+
#   geom_smooth()
# 
# ggplot(merge(cbind(gseaOrig$testData, preds=gseaOrig$testPreds), allAri), aes(x=value, y=as.numeric(gseaScaledImp)))+
#   geom_point()+
#   facet_wrap(~name)+
#   geom_smooth()
