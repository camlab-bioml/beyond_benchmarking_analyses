library(here)
library(ggplot2)

meta <- read.csv(here("collated_metadata.csv"), na.strings="")

num_unique <- lapply(meta, function(x){
  length(unique(x))
  }
  )
num_unique <- unlist(num_unique)

num_empty <- lapply(meta, function(x){
  length(which(is.na(x)))
})

num_empty <- as.data.frame(unlist(num_empty))

p1 <- ggplot(num_empty, aes(x=unlist(num_empty))) +
  geom_histogram(binwidth=5)+
  xlab("Number of missing entries")+
  cowplot::theme_cowplot()+
  ylab("Number of metadata fields")+
  ggtitle("Histogram of number of missing entries\nfor each metadata field")+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        title=element_text(size=40),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20))


num_unique <- as.data.frame(num_unique)
# p2 <- ggplot(num_unique, aes(x=unlist(num_unique))) +
#   geom_histogram(binwidth=5)+
#   xlab("Number of unique entries")+
#   cowplot::theme_cowplot()+
#   ylab("Number of metadata fields")+
#   ggtitle("Histogram of number of unique entries for each metadata field")+
#   theme(axis.title.x=element_text(size=30),
#         axis.title.y=element_text(size=30),
#         title=element_text(size=40))


# pdf(here("results", "figures", "supplementalFigs", "metadata_num_unique_hist.pdf"), height=15, width=15)
# print(p2)
# dev.off()

pdf(here("results", "figures", "supplementalFigs", "metadata_num_missing_hist.pdf"), height=15, width=15)
print(p1)
dev.off()


meta_stats <- cbind(num_empty, num_unique)
colnames(meta_stats) <- c("num_empty", "num_unique")
p3 <- ggplot(meta_stats, aes(x=num_empty, y=num_unique))+
  geom_point(size=4)+
  geom_jitter()+
  xlab("Number of missing entries")+
  ylab("Number of unique entries")+
  ggtitle("Number of unique entries vs. number of missing\nentries per metadata field")+
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        title=element_text(size=40),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20))

pdf(here("results", "figures", "supplementalFigs", "metadata_num_missing_unique_scatter.pdf"), height=16, width=16)
print(p3)
dev.off()

