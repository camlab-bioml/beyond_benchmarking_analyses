library(ggplot2)

singleR_ari <- read.csv("data/singleR/E-GEOD-114530-supervised-singleR.csv")

singleR_ari <- singleR_ari[,c("pipelines", "aris")]
model <- readRDS(here("data", "RFPreds_NumericParams.RDS"))
sil <- model$sil

data <- model$sil$testData
preds <- model$sil$testPreds

data <- cbind(data, preds)
data <- data[data$name == "E-GEOD-114530",]

data <- merge(data, singleR_ari, by="pipelines")

data <- data[order(-data$preds),]

data$rank <- 1:288


corARI <- round(cor(data$rank, data$aris), 3)

cor_df <- data.frame(aris=0.6, rank=1, cor=corARI)
silPlot <-ggplot(data, aes(x=aris, y=rank))+
  geom_point()+
  ylab("Predicted rank")+
  xlab("ARI with singleR annnotation")+
  cowplot::theme_cowplot()+
  theme(plot.title =
          element_text(hjust = 0.5, size=35),
        axis.title=element_text(size=15),
        axis.text=element_text(size=10))+
  ggtitle("SIL")



### DB ###
db <- model$db

data <- model$db$testData
preds <- model$db$testPreds

data <- cbind(data, preds)
data <- data[data$name == "E-GEOD-114530",]

data <- merge(data, singleR_ari, by="pipelines")

data <- data[order(-data$preds),]

data$rank <- 1:288


corARI <- round(cor(data$rank, data$aris), 3)

cor_df <- data.frame(aris=0.6, rank=1, cor=corARI)
dbPlot <-ggplot(data, aes(x=aris, y=rank))+
  geom_point()+
  ylab("Predicted rank")+
  xlab("ARI with singleR annnotation")+
  cowplot::theme_cowplot()+
  theme(plot.title =
          element_text(hjust = 0.5, size=35),
        axis.title=element_text(size=15),
        axis.text=element_text(size=10))+
  ggtitle("DB")

### CH ###
ch <- model$ch

data <- model$ch$testData
preds <- model$ch$testPreds

data <- cbind(data, preds)
data <- data[data$name == "E-GEOD-114530",]

data <- merge(data, singleR_ari, by="pipelines")

data <- data[order(-data$preds),]

data$rank <- 1:288


corARI <- round(cor(data$rank, data$aris), 3)

cor_df <- data.frame(aris=0.6, rank=1, cor=corARI)
chPlot <-ggplot(data, aes(x=aris, y=rank))+
  geom_point()+
  ylab("Predicted rank")+
  xlab("ARI with singleR annnotation")+
  cowplot::theme_cowplot()+
  theme(plot.title =
          element_text(hjust = 0.5, size=35),
        axis.title=element_text(size=15),
        axis.text=element_text(size=20))+
  ggtitle("CH")

### GSEA ###
gsea <- model$gsea

data <- model$gsea$testData
preds <- model$gsea$testPreds

data <- cbind(data, preds)
data <- data[data$name == "E-GEOD-114530",]

data <- merge(data, singleR_ari, by="pipelines")

data <- data[order(-data$preds),]

data$rank <- 1:288


corARI <- round(cor(data$rank, data$aris), 3)

cor_df <- data.frame(aris=0.6, rank=1, cor=corARI)
gseaPlot <-ggplot(data, aes(x=aris, y=rank))+
  geom_point()+
  ylab("Predicted rank")+
  xlab("ARI with singleR annnotation")+
  cowplot::theme_cowplot()+
  theme(plot.title =
          element_text(hjust = 0.5, size=35),
        axis.title=element_text(size=15),
        axis.text=element_text(size=10))+
  ggtitle("GSEA")


pdf("results/figures/singleR-comparison-ARI.pdf")
chPlot + dbPlot + silPlot + gseaPlot
dev.off()
