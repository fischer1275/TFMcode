# ---- Preliminaries ----
library("xlsx")
library("fastcluster")
library("gplots")
library("dendextend")
library(reshape)
library(matrixStats)
library("readxl")
library("heatmap3")
library("compare")
library("grid")
library("gridGraphics")
library("DESeq2")
library("ggplot2")
library("scater")
library("ggpubr")
library("GSVA")
library("RColorBrewer")
options(java.parameters = "-Xmx8g")

# ---- Update Annoations----

data=read.csv(file="MayoRNAseq_RNAseq_TCX_geneCounts_normalized.tsv", header=TRUE, sep="\t")

sampleLabels=substring(colnames(data)[-1],2)


convertSymbols=read_excel("symbolOut.xlsx")
ensID=convertSymbols$ensembl_gene_id
geneSymb=convertSymbols$hgnc_symbol

iA=match(ensID,as.character(data$ensembl_id))
iB=match(as.character(data$ensembl_id),ensID)

geneSymbSort=geneSymb[iB]
dataGeneSymb=data.frame(geneSymbSort,data)

 
covariates=read.csv(file="MayoRNAseq_RNAseq_TCX_covariates.csv", header=TRUE, sep=",")

iCov=match(sampleLabels,covariates$ID)

pheno=c("Diagnosis",NA,as.character(covariates$Diagnosis[iCov]))

dataOut=rbind((pheno),(dataGeneSymb))
dataOut=append(as.list(pheno),as.list(dataGeneSymb))

indCTRL=which(dataOut[1,]=="Control")
indAD=which(dataOut[1,]=="AD")

dataOut2=dataOut[,c(1:2,indCTRL,indAD)]
dataOut3=dataOut2[c(TRUE,!is.na(dataOut2[-1,1])),]


#iDuplicates=which(duplicated(dataOut3$geneSymbSort))

geneDuplicates=dataOut3$geneSymbSort[duplicated(dataOut3$geneSymbSort)]
uniqueDuplicates=unique(dataOut3$geneSymbSort[duplicated(dataOut3$geneSymbSort)])

iDuplicates= which(dataOut3$geneSymbSort %in% uniqueDuplicates)

indKeep=matrix(NA,nrow=1,ncol=length(uniqueDuplicates))
for (ind in 1:length(uniqueDuplicates)) 
{
  iDup=which(dataOut3$geneSymbSort==uniqueDuplicates[ind])
  dupData=dataOut3[iDup,3:dim(dataOut3)[2]]
  dupVar=apply(dupData,1,var)
  iVar=which.max(dupVar)
  
  indMaxVar=iDup[iVar]
  indKeep[ind]=indMaxVar
  
}

indRemove=setdiff(iDuplicates,indKeep)

dataOut4=dataOut3[-indRemove,]

write.csv(dataOut4,file="Mayo_Final.csv")



