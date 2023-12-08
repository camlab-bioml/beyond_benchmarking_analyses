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
  library(matrixStats)
  library(stringr)
  library(stringi)
  library(ggforce)
  library(umap)
})

plotUMAP <- function(train, test, metric){
  x <- rbind(train, test)
  x <- as_tibble(x)
  scores <- x %>%
    select(pipelines, name, ends_with("Imp"))%>%
    pivot_wider(names_from=name, values_from=ends_with("Imp")) %>%
    mutate(across(starts_with("E"), as.numeric))
  print(scores)
  UMAP <- umap(select(scores, starts_with("E")))
  
  scores <- scores %>%
    mutate(filt=case_when(grepl("lenient", pipelines) ~ "lenient",
                          grepl("default", pipelines) ~ "default",
                          grepl("stringent", pipelines) ~ "stringent")) %>%
    mutate(norm=case_when(grepl("scran", pipelines) ~ "scran",
                          grepl("sctransform", pipelines) ~ "sctransform",
                          grepl("seurat", pipelines) ~ "seurat"))%>%
    mutate(res=unlist(lapply(strsplit(pipelines, split="resolution."),"[",2)))%>%
    mutate(res=unlist(lapply(strsplit(res, split=".min.size"),"[",1)))%>%   
    mutate(dims=unlist(lapply(strsplit(pipelines, split="dims."),"[",2)))%>%
    mutate(dims=unlist(lapply(strsplit(dims, split=".k"),"[",1)))
  
  #print(scores$norm)
  UMAPdf <- as.data.frame(cbind(as.data.frame(UMAP$layout), filt=scores$filt, res=scores$res, norm=scores$norm, dims=scores$dims))
  
  return(UMAPdf)
}
designMatTest <- readRDS("corrected_RF_testMatNumeric_unscaled.RDS")
designMatTrain <- readRDS("corrected_RF_trainMatNumeric_unscaled.RDS")

silTrain <- cbind(designMatTrain$data, silCorImp = designMatTrain$labels$silCorImp)
silTest <- cbind(designMatTest$data, silCorImp =designMatTest$labels$silCorImp)

dbTrain <- cbind(designMatTrain$data, dbCorImp = designMatTrain$labels$dbCorImp)
dbTest <- cbind(designMatTest$data, dbCorImp =designMatTest$labels$dbCorImp)

chTrain <- cbind(designMatTrain$data, chCorImp = designMatTrain$labels$chCorImp)
chTest <- cbind(designMatTest$data, chCorImp =designMatTest$labels$chCorImp)

gseaTrain <- cbind(designMatTrain$data, gseaImp = designMatTrain$labels$gseaScaledImp)
gseaTest <- cbind(designMatTest$data, gseaImp =designMatTest$labels$gseaScaledImp)


silUMAP <- plotUMAP(silTrain, silTest, "Silhouette")
dbUMAP <- plotUMAP(dbTrain, dbTest, "DB")
chUMAP <- plotUMAP(chTrain, chTest, "CH")
gseaUMAP <- plotUMAP(gseaTrain, gseaTest, "GSEA")

UMAPdf <- bind_rows( CH=chUMAP, DB=dbUMAP, SIL=silUMAP, GSEA=gseaUMAP, .id="Metric")
UMAPdf$Metric <- factor(UMAPdf$Metric, levels=c("CH", "DB", "SIL", "GSEA"))


pdf("PaperFigures/pipelinesUMAP.pdf")
ggplot(UMAPdf, aes(x=as.numeric(V1), y=as.numeric(V2), color=as.factor(norm)))+
  geom_point()+
  xlab("UMAP1")+
  ylab("UMAP2")+
  cowplot::theme_cowplot()+
  #ggtitle(sprintf("%s UMAP",Metric))+ 
  facet_wrap(~Metric, scales="free")+
  scale_color_manual(values=c("scran" = "#CF597E", "sctransform" = "#EAE29C", "seurat" = "#089392" ))+
  guides(color=guide_legend(title="Normalization strategy"))+
  theme(strip.background = element_rect(fill="white"),
       strip.text = element_text(face="bold", size=20),
       axis.title.x=element_text(size=18),
       axis.title.y=element_text(size=18),
       legend.position = "bottom",
       legend.text = element_text(size=16))
dev.off()
