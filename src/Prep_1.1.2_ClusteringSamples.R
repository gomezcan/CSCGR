#############################################################################
# This  script make the kmeans clustering of the samples based on optimum   #
# number of clusters chosed, according to plots from Prep_1.1.1             #
#############################################################################

library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(factoextra)

# read temporar file from Prep_1.1.1
d_tsne_1_original <- read.table("Tem_tsne2D_DF.txt", h=T)

# set paratemers
nClusters= args[1]
nClusters= 12 # In this case work we used 12 clsuters

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1_original), nClusters) # 

# add cluster number to samples 
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Exploratory plot

# labels
clusterLabels <- d_tsne_1_original %>% 
  group_by(cl_kmeans) %>% select(V1, V2) %>% summarize_all(mean)

# Plot function
plot_cluster=function(data, var_cluster, palette, labels) {
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.35) +
    guides(colour=guide_legend(override.aes=list(size=3))) +
    xlab("tsne 1") + ylab("tsne 2") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_color_manual(values = rainbow(palette)) +
    geom_label_repel(data=labels, aes(label = cl_kmeans)) +
    guides(colour = FALSE) 
}

# Plot
plot_k <- plot_cluster(d_tsne_1_original, "cl_kmeans", nClusters, clusterLabels) 
# save plot
ggsave("t-SNE_plot_clusters.pdf")

## Save results
write.table(d_tsne_1_original, "tsne2D_clusters.txt", row.names = F, sep="\t", quote = F)

## remove tem file
system(paste("rm", "Tem_tsne2D_DF.txt", sep=' '))
