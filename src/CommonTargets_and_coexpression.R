library(scales)
library(ggrepel)
library(viridis)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggpubr)
library(plyr); library(dplyr)
library(data.table)

## Read data
# 1. Read results from common target results
# 2. Mapp co-expressed data into pairs. 

DataCommonTargets <- read_table2("CommonTarget_Resuslt.txt")
## Example file DataCommonTargets
# TF1	TF2	nTF1	nTF2	Common	J	Pval	OddR	Padj
# AT1G01060	AT1G01250	10374	403	151	0.014	0.292	1.0638	0.393
# AT1G01060	AT1G01720	10374	5650	2099	0.151	0.029	1.061	0.055
# AT1G01060	AT1G02230	10374	4055	1412	0.108	0.963	0.939	1
# AT1G01060	AT1G02250	10374	3784	1363	0.107	0.524	0.998	0.634


GetPCCdata2Net <- function(DFnet){
  TFs <- unique(DFnet$TF1)
  outdf <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 10)))
  colnames(outdf) <- c("TF1", "TF2", "nTF1", "nTF2", "Common","J","Pval","Padj","MR","PCC")
  c=1
  for (tf in TFs){
    File1 <- file.exists(paste("../DBcor_byGene/PCC.",tf,".txt", sep="")) # check if file in DB
    tem <- subset(DFnet, TF1==tf)  # subset TF1 pairs
    print(tf)
    if(File1=="TRUE"){
      TFs2 <- unique(tem$TF2)
      PCCdf1 <- as_tibble(fread(paste("../DBcor_byGene/PCC.",tf,".txt", sep=""), select = c(1:3), data.table = FALSE))
      PCCdf1 <- subset(PCCdf1, TAIRid %in% TFs2)
      colnames(PCCdf1)[c(1,2)] <- c("TF2","MR")
      tem <- left_join(tem, PCCdf1, by="TF2")
      outdf <- rbind(outdf, tem)
      print(paste("... TF:",tf,":",c,"Done...", sep=" "))
      c=c+1
    } else {
      print(paste("... TF:",tf,":",c,"No coexpress data...", sep=" "))
      tem[,"MR"] <- 0 # if not coexpression file available then MR and PCC equal to 0
      tem[,"PCC"] <- 0 # if not coexpression file available then MR and PCC equal to 0
      outdf <- rbind(outdf, tem)
      c=c+1
    }  
  }
  return(outdf)
}
# add PCC info
DataCommonTargets_PCC <- GetPCCdata2Net(DataCommonTargets)
# reduce data to threshols reported by Zand et al., 2018
DataCommonTargets_PCC_Results <- subset(DataCommonTargets_PCC, Padj <=0.00001 & PCC >=0.1)
write.table(DataCommonTargets_PCC_Results[,c(1,2,8,9,11)], "TFspairs_ResultsZand.txt", row.names = F, quote = F, sep="\t")
