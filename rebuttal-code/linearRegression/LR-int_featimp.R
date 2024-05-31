library(caret)

args <- commandArgs(trailingOnly=TRUE)
models <- readRDS(args[[1]])

sil <- models$sil
db <- models$db
ch <- models$ch
gsea <- models$gsea

allCoefs <- as.matrix(cbind(sil=sil$coefs, db=db$coefs, ch=ch$coefs, gsea=gsea$coefs))
colnames(allCoefs) <- c("SIL", "DB", "CH", "GSEA")

for (i in 1:length(rownames(allCoefs))){
  if (grepl("^V", rownames(allCoefs)[[i]])){
    num <- as.numeric(strsplit(rownames(allCoefs)[[i]], split="V")[[1]][[2]])-1
    rownames(allCoefs)[[i]] <- paste0("PC", num)
  }
}

allCoefs <- allCoefs[-which(rowSums(allCoefs)==0),]
allCoefs <- allCoefs[-which(rownames(allCoefs)=="(Intercept)"),]
print(allCoefs)
n=nrow(allCoefs)*ncol(allCoefs)
f2 = circlize::colorRamp2(seq(-0.01, 0.01, length = nrow(allCoefs)), hcl.colors(nrow(allCoefs),"Purple-Green"))

# pdf("PaperFigures/LRCoefficientsHeatmap.pdf")
pdf("/home/campbell/cfang/bb-rebuttal/results/figures/linearRegression/LR-int_featImp.pdf", height=30, width=20)
ComplexHeatmap::Heatmap(allCoefs, col=f2, name="Feature\ncoefficient",
                        row_title = "Feature",
                        column_title="Metric",
                        column_title_side="bottom")
dev.off()
