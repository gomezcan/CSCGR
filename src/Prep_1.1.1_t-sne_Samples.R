########################################################################################### 
# This script recive as input the expression matrix to calculated the t-sne values. Then, #
# using the tsne values make and evaluatoin of an optimum of cluster determinen by the    #
# within cluster sums of squares. 
###########################################################################################

library(RColorBrewer)
library(Rtsne)
library(factoextra)


set.seed(1234)  

# Read Expression matrix
ExpressionDB <- read.table("Ath_r_v15_08_expression/Ath_r_v15_08_expressionTAIRids.txt.gz", h=T)
# GeneID as row names
row.names(ExpressionDB) <- ExpressionDB$TAIRid
# Remove GeneID and ENTREZ ID
ExpressionDB <- ExpressionDB[,-c(1,1417)]


##################   Calculate t-sne values    ##########################

tsne_model_1 = Rtsne(as.matrix(t(ExpressionDB)), check_duplicates=FALSE, 
                     pca=TRUE, perplexity=30, theta=0.5, dims=2)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y) 


## Exploratory plot 1

Plot_1 <- ggplot(d_tsne_1, aes(x=V1, y=V2)) +  
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6)))+
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")
ggsave("t-SNE_plot.pdf")

#### Test number of kmeans
## Exploratory plot 2 
fviz_nbclust(d_tsne_1, kmeans, method = "wss", nboot = 300, k.max = 25)
ggsave("Optimum_bumber_of_clusters.pdf")

## save tem file
d_tsne_1$Sample <- names(ExpressionDB) # add samples ID
write.table(d_tsne_1, "Tem_tsne2D_DF.txt", row.names = F, sep="\t", quote = F)

