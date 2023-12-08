suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(ggplot2)
  library(Metrics)
  library(gridExtra)
  library(caret)
  library(matrixStats)
  library(randomForest)
  #library(pcaMethods)
  library(stringr)
  library(stringi)
  library(glmnet)
  library(ggforce)
})
# Plot the metrics vs. number of clusters
designMatCor <- readRDS("corrected_designMat_unscaled.RDS")
designMatUncorr <- readRDS("uncorrected_designMatrix_scaled_cleaned.RDS")
ari <- readRDS("ari_unscaled_cleaned.RDS")
dataset = "E-HCAD-25"

datasetCor <- designMatCor %>% 
  filter(.,grepl(dataset, name)) %>%
  select(nclusts, silCorImp, dbCorImp, chCorImp, filt, norm, dims) %>%
  rename(sil=silCorImp) %>%
  rename(db=dbCorImp) %>%
  rename(ch=chCorImp) %>%
  mutate(across(starts_with(c("sil", "db", "ch", "nclusts", "dims")),as.numeric))

datasetUncorr <- designMatUncorr %>%
  filter(., grepl(dataset, name))%>%
  select(nclusts, sil, db, ch, filt, norm, dims)%>%
  mutate(db=db*-1) %>%
  mutate(across(starts_with("dims"), as.numeric))# because corrected db had sign flipped before correction

df <- bind_rows(list("Uncorrected"=datasetUncorr, "Corrected"=datasetCor), .id="cor")
df <- df %>%
  pivot_longer(cols=c(sil, db,ch), names_to="metric")%>%
  mutate(metric=toupper(metric))

df$cor <- factor(df$cor, levels=c("Uncorrected", "Corrected"))

pdf(sprintf("PaperFigures/%sMetricsNclusts.pdf", dataset))
ggplot(df, aes(x=as.numeric(nclusts), y=as.numeric(value)))+
  geom_point()+
  ylab("Metric value on dataset E-HCAD-25")+
  xlab("Number of clusters found by pipeline on dataset E-HCAD-25")+
  facet_grid(metric~cor, scales="free")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 1)+
  geom_smooth()
dev.off()

pdf(sprintf("PaperFigures/%sMetricsNclustsFilt.pdf", dataset))
ggplot(df, aes(x=filt, y=as.numeric(value)))+
  geom_point(alpha=0.15)+
  ylab("Metric value on dataset E-HCAD-25")+
  xlab("Filtering strategy applied to dataset E-HCAD-25")+
  facet_grid(metric~cor, scales="free")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 1)
dev.off()

pdf(sprintf("PaperFigures/%sMetricsNclustsNorm.pdf", dataset))
ggplot(df, aes(x=norm, y=as.numeric(value)))+
  geom_point(alpha=0.15)+
  ylab("Metric value on dataset E-HCAD-25")+
  xlab("Normalization strategy applied to dataset E-HCAD-25")+
  facet_grid(metric~cor, scales="free")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 1)
dev.off()

pdf(sprintf("PaperFigures/%sMetricsNclustsDims.pdf", dataset))
ggplot(df, aes(x=as.numeric(dims), y=as.numeric(value)))+
  geom_point(alpha=0.15)+
  ylab("Metric value on dataset E-HCAD-25")+
  xlab("# of dimensions for dimensionality reduction on dataset E-HCAD-25")+
  facet_grid(metric~cor, scales="free")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 1)
dev.off()

gseaDf <- designMatCor %>% 
  select(nclusts, gseaScaledImp, name, pipelines) %>%
  rename(GSEA=gseaScaledImp)

# pdf("PaperFigures/gseaNclusts.pdf")
# for (i in 1:8){
#   p <- ggplot(gseaDf, aes(x=as.numeric(nclusts), y=as.numeric(GSEA)))+
#     geom_point()+
#     ylab("Mean normalized enrichment score")+
#     xlab("Number of clusters found by pipeline")+
#     cowplot::theme_cowplot()+
#     theme(strip.background = element_rect(fill="white"),
#           strip.text = element_text(face="bold"))+
#     theme(plot.title = element_text(hjust = 0.5))+
#     facet_wrap_paginate(~name, scales="free", page=i, nrow=4, ncol=3)+
#     geom_smooth()
#   print(p)
# }
# dev.off()

allMetricsdf <- designMatUncorr %>%
  select(nclusts, sil, db, ch, name, pipelines, filt, norm, dims)%>%
  mutate(db=db*-1) # because corrected db had sign flipped before correction

allMetricsdf <- merge(allMetricsdf, gseaDf, by=c("name", "pipelines"))
allMetricsdf$GSEA <- as.numeric(allMetricsdf$GSEA)
allMetricsdf <- pivot_longer(allMetricsdf, cols = c(sil, db, ch, GSEA), names_to="metric")
allMetricsdf$metric <- toupper(allMetricsdf$metric)

allMetricsdf$metric_f <- factor(allMetricsdf$metric, levels=c("CH", "DB", "SIL", "GSEA"))

pdf("PaperFigures/uncorrectedMetricsNclusts.pdf")
ggplot(allMetricsdf, aes(x=as.numeric(nclusts.x), y=as.numeric(value)))+
  geom_point(alpha=0.15)+
  ylab("Metric value")+
  xlab("Number of clusters found by pipeline")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~ metric_f, scales="free")+
  geom_smooth()
dev.off()

pdf("PaperFigures/uncorrectedMetricsFilt.pdf")
ggplot(allMetricsdf, aes(x=filt, y=as.numeric(value)))+
  geom_point(alpha=0.15)+
  ylab("Metric value")+
  xlab("Filtering strategy")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~ metric_f, scales="free")
dev.off()

pdf("PaperFigures/uncorrectedMetricsNorm.pdf")
ggplot(allMetricsdf, aes(x=norm, y=as.numeric(value)))+
  geom_point(alpha=0.15)+
  ylab("Metric value")+
  xlab("Normalization strategy")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~ metric_f, scales="free")
dev.off()

pdf("PaperFigures/uncorrectedMetricsDims.pdf")
ggplot(allMetricsdf, aes(x=as.numeric(dims), y=as.numeric(value)))+
  geom_point(alpha=0.15)+
  ylab("Metric value")+
  xlab("Number of dimensions")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~ metric_f, scales="free")
dev.off()

sildf <- allMetricsdf %>%
  filter(metric=="SIL") %>%
  filter(name %in% c("E-CURD-10", "E-GEOD-89232","E-GEOD-75367","E-GEOD-83139","E-GEOD-125970", "E-ENAD-27"))

ggplot(sildf, aes(x=as.numeric(nclusts.x), y=as.numeric(value)))+
  geom_point(alpha=0.75)+
  ylab("Metric value")+
  xlab("Number of clusters found by pipeline")+
  cowplot::theme_cowplot()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~ name, scales="free")+
  geom_smooth()

ari_tbl <- ari %>%
  as_tibble() %>%
  mutate(across(starts_with("E"), scale)) %>%
  pivot_longer(cols=starts_with("E"), names_to="name", values_to="ari")
plotDf <- merge(designMatCor, ari_tbl, by=c("pipelines", "name"))

ggplot(plotDf, aes(x=ari, y=as.numeric(silCorImp)))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~name)

nclust <- designMatUncorr %>%
  select(pipelines, name, nclusts)
plotDf <- merge(plotDf, nclust, by=c("pipelines", "name"))

ariMat <- scale(as.matrix(ari[,-1]))
n=nrow(ariMat)*ncol(ariMat)
colFun = circlize::colorRamp2(seq(-4, 2, length = n), hcl.colors(n,"viridis"))
ariHM <- ComplexHeatmap::Heatmap(ariMat, col=colFun, 
                                 column_names_gp = grid::gpar(fontsize = 0),
                                 name="ARI")
pdf("PaperFigures/ariHM.pdf", width=5, height=7)
ComplexHeatmap::draw(ariHM, 
     row_title="Pipeline",
     column_title_side="bottom",
     column_title="Dataset")
dev.off()

pdf("PaperFigures/ariNclusts.pdf")
for (i in 1:3){
  p<- ggplot(plotDf, aes(x=nclusts.x, y=plotDf$ari))+
    geom_point()+
    geom_smooth()+
    cowplot::theme_cowplot()+
    facet_wrap_paginate(~name, scales="free", page=i, nrow=3, ncol=2)+
    theme(strip.background = element_rect(fill="white"),
          strip.text = element_text(face="bold"))+
    ylab("ARI")+
    xlab("Number of clusters found by pipeline")
  print(p)
}
dev.off()

plotDfCor <- plotDf %>%
  group_by(name)%>%
  group_modify(~ as.data.frame(cbind(silAriCor = cor(.$ari, as.numeric(.$silCorImp)))))


pdf("PaperFigures/rfTuning.pdf")
plot(ch$model, main="CH")
plot(db$model, main="DB")
plot(gsea$model, main="GSEA")
plot(sil$model, main="SIL")
dev.off()
