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

data=read_excel("Zhang2013Data.xlsx", sheet="GSE44770_series_matrix", range="A51:HW39383")
sampleLabels=colnames(data)[-1]

dataOut=data[-seq(1:52),]
colnames(dataOut)=c(as.character(data[52,]))

reddphenoSource=data[20,-1]
pheno=phenoSource
pheno[phenoSource=="disease: A"]="AD"
pheno[phenoSource=="disease: N"]="CTRL"



rosettaAnnotations=read_excel("RosettaProbe.xlsx", sheet="GPL4372", range="A28:C39330")

iB=match(dataOut$ID_REF,rosettaAnnotations$ID)

geneSymbSort=rosettaAnnotations$`Gene symbol`[iB]
dataGeneSymb=data.frame(geneSymbSort,dataOut)

dataOut1=rbind(c(NA,NA,t(pheno)),dataGeneSymb)

indCTRL=which(dataOut1[1,]=="CTRL")
indAD=which(dataOut1[1,]=="AD")

dataOut2=dataOut1[,c(1:2,indCTRL,indAD)]
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

write.csv(dataOut4,file="MtSinai_Final.csv")
 


