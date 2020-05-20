#############################################################################
# This script calculate wPCC for a subset of samples defined based on the   #
# clustering analysis in Prep_1.1      					    #
# With  the dataset used in this work, this number was set to 12 by default #
#############################################################################

library(ggplot2)
library(reshape)
library(wCorr)
library(parallel)

###### Read input files ######
## Expresion data
ExpressionDB <- read.table("Ath_r_v15_08_expressionTAIRids.txt.gz", h=T)
# GeneID as row names
row.names(ExpressionDB) <- ExpressionDB$TAIRid
# Remove GeneID and ENTREZ ID
ExpressionDB <- ExpressionDB[,-c(1,1417)]

## Read Sample cluster  
tsne2d <- read.table("tsne2D_clusters.txt", h=T)

###### create directory to save results ########
system(paste("mkdir", "Local_CoExpDB", sep=' '))
###############################################

##

#################### Funtions ######################

# samples weight

SampleW <- function(expDFmod) {
 # Calculate the correlation matrix between samples within each condition
 # if the PCC is smaller than 0.4, set w to 0.4.
 # Note: 0.4 is the optimal threshold using by ATTED-II
  CexpDFmod <- as.data.frame(cor(expDFmod))
  CexpDFmod[CexpDFmod<0.4] <- 0.4 # if cor small than threshold then set to threshold
  CexpDFmod$redundancy <- apply((CexpDFmod-0.4)/(1-0.4),1,sum) 
  CexpDFmod$weight <- 1/sqrt(CexpDFmod$redundancy)
  return(CexpDFmod$weight)
}


# calculate the redundancy score and weight for samples within condition 
wCorRows <- function(id){
 # This function calcuate the wPCC of a Gene define by "id" vs all genes in Expressoin matrix
 # mdf: expression matrix of a sibset of samples
  wPCC_vector <- vector()
  for (i in 1:nrow(mdf)){
    wPCC_vector[i] <- weightedCorr(x=mdf[id,], y=mdf[i,], wvector, method = 'Pearson')
  }
  return(wPCC_vector) # vector of  wPCC of gene "id" vs all genes
}

#####################################################

################## Coexpression by Module data ################

Modules <-  unique(tsne2d$cl_kmeans)


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

  ## Calculated weight vector
  wvector <- SampleW(mdf)

  ## calcule wPCC parallelizing number of genes been analyzed
  Cm <- mclapply(row.names(mdf), wCorRows, mc.cores=10)

  ## reshape like square matrix 
  Cm <- data.frame(matrix(unlist(Cm), nrow=nrow(mdf), byrow=T),stringsAsFactors=FALSE)
  # Cm <- Cm[1:nrow(Cm), 1:ncol(Cm)]
  ## Clean matrix s
  Cm[is.na(Cm)] <- 0
  Cm[Cm == 1] <- NA
  Cm[is.na(Cm)] <- 0
  diag(Cm) <- 1
  # add gene ids
  colnames(Cm) <- row.names(mdf) 
  row.names(Cm) <- row.names(mdf) 
  ## print(paste(".......Donde Caculation : ",name, " ....", sep=''))
  
  # Rank values: large wPCC values (good) with "-" will be rank close to 1. 
  Rankm <- apply(-Cm, 1, rank)

  # Mutual rank calcuation 
  MRm <- round(sqrt(Rankm*t(Rankm)),2)

  # list of genes and files to write
  GeneID <- as.character(row.names(mdf))
  
  #print(paste(".......Saving: ",name, " ....", sep=''))
  for (id in GeneID){
    # create a DF by genes, and save it into Local_CoExpDB
    TemDF <- as.data.frame(cbind(PCC=Cm[,id], MR=MRm[,id]))
    TemDF <- TemDF[order(TemDF$MR),]
    TemDF["GeneID"] <- row.names(TemDF)
    write.table(TemDF, paste("Local_CoExpDB/PCC.",id,".",name,".","txt",sep=""), row.names = F, sep='\t', quote = F)
    
  }
 
rm(Cm)
rm(Rankm)
rm(MRm)
