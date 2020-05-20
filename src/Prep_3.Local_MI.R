#############################################################################
# This script calculate MI for a subset of samples defined based on the     #
# clustering analysis in Prep_1.1      					    #
# With  the dataset used in this work, this number was set to 12 by default #
#############################################################################

library(parmigene)

###### Read input files ######
## Expresion data
ExpressionDB <- read.table("../Ath_r_v15_08_expressionTAIRids.txt.gz", h=T)
# GeneID as row names
row.names(ExpressionDB) <- ExpressionDB$TAIRid
# Remove GeneID and ENTREZ ID
ExpressionDB <- ExpressionDB[,-c(1,1417)]

## Read Sample cluster  
tsne2d <- read.table("tsne2D_clusters.txt", h=T)

###### create directory to save results ########
system(paste("mkdir", "Local_CoExpDB_MI", sep=' '))
###############################################


################## Coexpression by Module data ################

Modules <-  unique(tsne2d$cl_kmeans)


## define THREADS to parallelized calculations 
OMP_NUM_THREADS=25



for (i in 1:Modules){
	# loop for walk over modules, make subset of expression data,
	# calculate wPCC, calculate mutual rank (MR), and save a file
	# each gene indicating the cluster number
	
	## Samples in clusters	
	S=as.character(subset(tsne2d, cl_kmeans==i)$Sample)
	name=paste("Module_", i,sep="")
	
	## Subsample expression matrix
	mdf <- as.matrix((ExpressionDB[,S]))
	print(paste(".......Caculating : ",name, " ....", sep=''))

	# Calculate MI
	MI <- knnmi.all(as.matrix(mdf), k=3, noise=1e-12)

	# Rank from smaller to bigger: large MI values (good) with "-" will be rank close to 1. 
	MI_Rank <- apply(-MI, 1, rank)

	# Mutual rank calcuation   
	MR_MI <- sqrt(MI_Rank*t(MI_Rank))
	MR_MI <- round(MR_MI, 2)

	GeneID <- as.character(row.names(ExpressionDB))

	## Save data ##
	for (i in GeneID){
		TemDF <- as.data.frame(cbind(MI=MI_Rank[,i], MR=MR_MI[,i]))
		TemDF <- TemDF[order(TemDF$MI),]
		TemDF <- cbind(GeneID=row.names(TemDF), TemDF)
		write.table(TemDF, paste("Local_CoExpDB_MI/MI.",id,".",name,".","txt",sep=""), row.names = F, sep='\t', quote = F)
}
}
