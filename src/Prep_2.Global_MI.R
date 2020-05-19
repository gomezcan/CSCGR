############################################################################
# This script recive expression data to calculate co-expresion based on MI #
############################################################################

library(parmigene)
library(dplyr)

###### Read input files ######
ExpressionDB <- read.table("../Ath_r_v15_08_expressionTAIRids.txt", h=T)
# GeneID as row names
row.names(ExpressionDB) <- ExpressionDB$TAIRid
# Remove GeneID and ENTREZ ID
ExpressionDB <- ExpressionDB[,-c(1,1417)]

### create directory to save results ###
system(paste("mkdir", "MI_CoExpDB", sep=' '))
########################################

## define THREADS to parallelized calculations 
OMP_NUM_THREADS=25

# Calculate MI
MI <- knnmi.all(as.matrix(ExpressionDB), k=3, noise=1e-12)

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
  write.table(TemDF, paste("MI_CoExpDB/MI.",i,".txt",sep=""), row.names = F, sep='\t', quote = F)
}


