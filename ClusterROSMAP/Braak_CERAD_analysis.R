options(java.parameters = "-Xmx8g")
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
#library("scater")
library("ggpubr")
library("GSVA")
library("RColorBrewer")
library("impute")
library("limma")
library("vsn")
library("Hmisc")
library("org.Hs.eg.db")
library("clusterProfiler")
library("qusage")
library("ReactomePA")
library("msigdbr")
library(AnnotationDbi)
source("../Custom R Functions/Mode.R")
source("../Custom R Functions/summarySE.R")
source("../Custom R Functions/rcorr.adjustFDR.R")



#addition 

library(factoextra)
library(cluster)
library("enrichplot")
library("rrvgo")

#Set heatmap3 color bar parameters
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red")

breakBarColorsCor=c(-200,seq(-1, 1, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColorsCor = colorpanel(length(breakBarColorsCor)-1, "black", "white", "orange")

# ---- Read Stanford Cell Type Data and Clustergram----
cellTypeDataCTD=read_excel("../Zhang2016GeneExpression_ForClustergram.xlsx", sheet="ClusterGeneSortThresh")
yCTD=cellTypeDataCTD[1:11077,2:18];
geneLabelsCTD=cellTypeDataCTD$Gene;


yMatCTD=matrix(unlist(yCTD),ncol=17,byrow=FALSE)

yMatCTDZ=t(apply(yMatCTD,1,scale)); #z-score the data


hrCTD= hclust(dist((yMatCTDZ),method = "euclidean"), method="average")
hcCTD= hclust(dist(t(yMatCTDZ),method = "euclidean"), method="average")


cutThreshCTD=3.76 #use with the euclidian distance method, average

myclCTD = cutree(hrCTD, h=cutThreshCTD); 
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Spectral"))
mycolhcCTD =  sample(colorBrewerInterp(length(unique(myclCTD))))
mycolhcCTD = mycolhcCTD[as.vector(myclCTD)]

t1CTD=as.dendrogram(hrCTD)

#par(oma=c(5,15,1,1))
heatmap3(yMatCTDZ,  
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=as.dendrogram(hrCTD), Colv=as.dendrogram(hcCTD),  scale="none",
         RowSideColors=mycolhcCTD, RowSideLabs=NA, cexRow=3,cexCol =1,
         labRow = NA, labCol = colnames(yCTD), 
         margins = c(5,20), ColSideWidth=.4)


rowOrder=order.dendrogram(t1CTD);
rowOrderCol=order.dendrogram(as.dendrogram(hcCTD));


#Prepare clustered data for writing to table
yMatClust=yMatCTD[rowOrder,rowOrderCol];
geneLabelsClustCTD=geneLabelsCTD[rowOrder];
clusterNumberCTD=myclCTD[rowOrder];


# Plot the astrocyte cluster 4


heatmap3(yMatClust[clusterNumberCTD==4,],  
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = colnames(yCTD)[rowOrderCol])

# Plot the neuron cluster 13


heatmap3(yMatClust[clusterNumberCTD==13,],  
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = colnames(yCTD)[rowOrderCol])


#plot the microglia cluster 17 

heatmap3(yMatClust[clusterNumberCTD==17,],  
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = colnames(yCTD)[rowOrderCol])


# ---- Read and Plot Mt Sinai Data ----

mtSinaiData=read.csv("./Data Pre-Processing/ROSMAP_Final.csv", sep=",")

pathScoresData=read.csv("./Data Pre-Processing/ROSMAP_pathScores.csv", sep=",")


y=mtSinaiData[2:dim(mtSinaiData)[1],4:dim(mtSinaiData)[2]]; #Careful of indexing
geneLabelsMtSinai=mtSinaiData$geneSymbSort[2:nrow(mtSinaiData)];    #Careful of indexing 
phenoMtSinai=(mtSinaiData[1,4:dim(mtSinaiData)[2]])
pathScores=pathScoresData[,4:ncol(pathScoresData)]

#Re-write CERAD into standard scale
CERADTemp=pathScores[2,]
pathScores[2,CERADTemp==4]=0
pathScores[2,CERADTemp==3]=1
pathScores[2,CERADTemp==2]=2
pathScores[2,CERADTemp==1]=3

rownames(pathScores)=c("braak","cerad","apoe","sex","pmi")


phenoMtSinai <- data.frame(lapply(phenoMtSinai, as.character), stringsAsFactors=FALSE)

yMat=apply(matrix(unlist(y),ncol=dim(y)[2],byrow=FALSE),2,as.numeric)

#Impute missing values
yMatImp=impute.knn(yMat,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
yMat=yMatImp$data;
colnames(yMat)=colnames(phenoMtSinai)


colorBar=ifelse(phenoMtSinai=="CTRL","navy","red3")
colorBar[phenoMtSinai=="MCI"]="peachpuff2";


heatmap3(yMat, ColSideColors =c(colorBar), ColSideWidth=1, 
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA) 


# ---- Select the genes in AD dataset according to CTD and cluster, plot----
iMatch=match(geneLabelsClustCTD,geneLabelsMtSinai)

#ySort=matrix(data = NA, nrow = dim(yMatCTD)[1], ncol = dim(yMat)[2]);
ySort=yMat[iMatch,];
geneSort=geneLabelsMtSinai[iMatch]

#Plotting sorted genes
#Put everything in same heatmap
ySortZ=t(apply(ySort,1,scale)); #z-score the data
colnames(ySortZ)=colnames(ySort)
#ySortZNeurons[which(is.na(ySortZNeurons))]=0; #replaces NAs due to constant columns with zeros.
#No clustering

#Column clustering whole dataset
hcSort= hclust(dist(t(ySortZ),method = "euclidean"), method="ward.D2")
heatmap3(ySortZ, ColSideColors =c(colorBar), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSort),  scale="none",
         labRow = NA, labCol = NA) 

#Separate heatmaps and clustering
ySortZCTRL=ySortZ[,phenoMtSinai=="CTRL"];
ySortZMCI=ySortZ[,phenoMtSinai=="MCI"];
ySortZAD=ySortZ[,phenoMtSinai=="AD"];


hcSortCTRL= hclust(dist(t(ySortZCTRL),method = "euclidean"), method="ward.D2")
heatmap3(ySortZCTRL, ColSideColors =c(colorBar)[1:sum(phenoMtSinai=="CTRL")], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortCTRL),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 

hcSortMCI = hclust(dist(t(ySortZMCI),method = "euclidean"), method="ward.D2")
heatmap3(ySortZMCI, ColSideColors =c(colorBar)[sum(phenoMtSinai=="CTRL")+1:sum(phenoMtSinai!="AD")], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortMCI),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 

hcSortAD= hclust(dist(t(ySortZAD),method = "euclidean"), method="ward.D2")
heatmap3(ySortZAD, ColSideColors =c(colorBar)[sum(phenoMtSinai!="AD")+1:length(phenoMtSinai)], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortAD),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 




#----Cell Type Normalization----
#Setup phenotype used for DEseq2 format below
phenoNum=as.integer(ifelse(phenoMtSinai=="CTRL",0,1))
phenoNum[phenoMtSinai=="AD"]=2
condition=as.factor(c(phenoNum))

#Clean up data - keep rows that are not NA
indGenesKeep=which(!is.na(ySort[,1]))
ySortClean=ySort[indGenesKeep,]
geneSortClean=geneSort[indGenesKeep]
clusterNumberCTDClean=clusterNumberCTD[indGenesKeep]  

#Initialize normalized data
ySortNorm=ySortClean

#create matrix with gene annotation for later extraction of genes in various subgroups
ySortNorm_ann <- ySortClean 
rownames(ySortNorm_ann) <- geneSortClean


#Normalize neurons
ySortNormNeuronGenes=normalizeVSN(ySortNorm[clusterNumberCTDClean==13,])
ySortNeuronGenes=ySortClean[clusterNumberCTDClean==13,]
ySortNorm[clusterNumberCTDClean==13,]=normalizeVSN(ySortNorm[clusterNumberCTDClean==13,])

neuroGeneNames <- rownames(ySortNorm_ann[clusterNumberCTDClean==13,])


#Normalize astrocytes
ySortNormAstrocyteGenes=normalizeVSN(ySortNorm[clusterNumberCTDClean==4,])
ySortAstrocyteGenes=ySortNorm[clusterNumberCTDClean==4,]
ySortNorm[clusterNumberCTDClean==4,]=normalizeVSN(ySortNorm[clusterNumberCTDClean==4,])

astroGeneNames <- rownames(ySortNorm_ann[clusterNumberCTDClean==4,])

#Normalize microglia
ySortNormMicrogliaGenes=normalizeVSN(ySortNorm[clusterNumberCTDClean==17,])
ySortMicrogliaGenes=ySortNorm[clusterNumberCTDClean==17,]
ySortNorm[clusterNumberCTDClean==17,]=normalizeVSN(ySortNorm[clusterNumberCTDClean==17,])

microGeneNames <- rownames(ySortNorm_ann[clusterNumberCTDClean==17,])

#Z score neuron gene normalized for plotting below
ySortNormZ=t(apply(ySortNorm,1,scale)); #z-score the data
ySortNormZCTRL=ySortNormZ[,phenoMtSinai=="CTRL"];
ySortNormZMCI=ySortNormZ[,phenoMtSinai=="MCI"];
ySortNormZAD=ySortNormZ[,phenoMtSinai=="AD"];



#Select sample indices from CTRL and AD Clustering Above (non-normalized)

#CTRL

hcCTRLCut = cutree(hcSortCTRL, k=2); 

colorhcCTRL = rainbow(length(unique(hcCTRLCut)), start=0.1, end=0.9); 
colorhcCTRL=c("skyblue4","skyblue")
colorhcCTRL = colorhcCTRL[as.vector(hcCTRLCut)]

startColorIndCTRL=1
endColorIndCTRL=sum(phenoMtSinai=="CTRL")

combinedColorCTRL=cbind(colorBar[startColorIndCTRL:endColorIndCTRL], 
                        colorhcCTRL)
heatmap3(ySortZCTRL, ColSideColors = combinedColorCTRL, ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortCTRL),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 

heatmap3(ySortNormZCTRL, ColSideColors = combinedColorCTRL, ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortCTRL),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 



colOrderCTRL=order.dendrogram(as.dendrogram(hcSortCTRL));
clusterNumerCTRLSort=hcCTRLCut[colOrderCTRL] #Shows we want cluster 1

indCTRL_G1=which(hcCTRLCut==1)




#MCI

hcMCICut = cutree(hcSortMCI, k=2); 

colorhcMCI=c("papayawhip","peachpuff3")
colorhcMCI = colorhcMCI[as.vector(hcMCICut)]

startColorIndMCI=sum(phenoMtSinai=="CTRL")+1
endColorIndMCI=sum(phenoMtSinai=="CTRL")+sum(phenoMtSinai=="MCI")

combinedColorMCI=cbind(colorBar[startColorIndMCI:endColorIndMCI], 
                       colorhcMCI)
heatmap3(ySortZMCI, ColSideColors = combinedColorMCI, ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortMCI),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 


colOrderMCI=order.dendrogram(as.dendrogram(hcSortMCI));
clusterNumerMCISort=hcMCICut[colOrderMCI] #Shows we want cluster 1

indMCI_G1=which(hcMCICut==1)


#AD
hcADCut = cutree(hcSortAD, k=2); 
colorhcAD = rainbow(length(unique(hcADCut)), start=0.3, end=0.6); 
colorhcAD=c("orangered","orangered4")
colorhcAD = colorhcAD[as.vector(hcADCut)]

startColorIndAD=sum(phenoMtSinai!="AD")+1
endColorIndAD=length(phenoMtSinai)

combinedColorAD=cbind(colorBar[startColorIndAD:endColorIndAD], 
                      colorhcAD)
heatmap3(ySortZAD, ColSideColors = combinedColorAD, ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortAD),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 

heatmap3(ySortNormZAD, ColSideColors = combinedColorAD, ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortAD),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 


colOrderAD=order.dendrogram(as.dendrogram(hcSortAD));
clusterNumerADSort=hcADCut[colOrderAD] #Shows we want cluster 1

indAD_G1=which(hcADCut==1)


#Define Colors for each group
combinedColorALL=cbind(c(colorBar[startColorIndCTRL:endColorIndCTRL], colorBar[startColorIndMCI:endColorIndMCI], colorBar[startColorIndAD:endColorIndAD]), 
                       c(colorhcCTRL,colorhcMCI, colorhcAD))

allColorGroups=c(colorhcCTRL,colorhcMCI,colorhcAD)

indBlue=which(c(colorhcCTRL,colorhcMCI,colorhcAD)=="skyblue4") #colors for the 2 control subgroups 
indGreen=which(c(colorhcCTRL,colorhcMCI,colorhcAD)=="skyblue")
indPlum=which(c(colorhcCTRL,colorhcMCI,colorhcAD)=="papayawhip") #colors for the 2 MCI subgroups 
indIndianred=which(c(colorhcCTRL,colorhcMCI,colorhcAD)=="peachpuff3")
indRed=which(c(colorhcCTRL,colorhcMCI,colorhcAD)=="orangered") #colors for the 2 AD subgroups 
indPurple=which(c(colorhcCTRL,colorhcMCI,colorhcAD)=="orangered4")

indSort=c(indBlue,indGreen,indPlum,indIndianred,indRed,indPurple)
allColorGroups=allColorGroups[indSort]
#geneSetEnrichSort=geneSetEnrich[,indSort]


#Write out neuron genes sorted with samples sorted by indSort
yNeuROSMAP=ySortNormNeuronGenes[,indSort]
rownames(yneuROSMAP)=geneSortClean[clusterNumberCTDClean==13]


#Write out astrocyte genes sorted with samples sorted by indSort
yAstROSMAP=ySortNormAstrocyteGenes[,indSort]
rownames(yAstROSMAP)=geneSortClean[clusterNumberCTDClean==4]



#Write out astrocyte genes sorted with samples sorted by indSort
yMicROSMAP=ySortNormMicrogliaGenes[,indSort]
rownames(yMicROSMAP)=geneSortClean[clusterNumberCTDClean==17]

phenoROSMAP=allColorGroups
#save(yAstROSMAP, phenoROSMAP, file="../ADGeneSignature/SigROSMAP.RData")


#----FC----
#Compute histogram comparing fold change between normlized and non-normlized genes

#Picking groups 1 for CTRL and AD
indCTRL_FC=which(phenoMtSinai=="CTRL" & c(hcCTRLCut==1,rep(NA,endColorIndAD-startColorIndMCI+1)))
indAD_FC=which(phenoMtSinai=="AD" & c(rep(NA,endColorIndMCI), hcADCut==1))


#Compute histogram based on normalized genes
#Neurons

ySortNeuronGenesZ=t(apply(ySortNeuronGenes,1,scale)); #z-score the data
ySortNormNeuronGenesZ=t(apply(ySortNormNeuronGenes,1,scale)); #z-score the data

meanCTRL=apply(ySortNeuronGenesZ[,indBlue],1,mean)
meanAD=apply(ySortNeuronGenesZ[,indPurple],1,mean)
logFC1=((meanAD)-(meanCTRL))

meanCTRL=apply(ySortNormNeuronGenesZ[,indBlue],1,mean)
meanAD=apply(ySortNormNeuronGenesZ[,indBlue],1,mean)
logFC2=((meanAD)-(meanCTRL))

neuoronLogFC=data.frame(
  logFC=c(logFC1,logFC2),
  norm=rep(c("Not Normalized","Normalized"),each=length(logFC1)))


#Astrocytes
ySortAstrocyteGenesZ=t(apply(ySortAstrocyteGenes,1,scale)); #z-score the data
ySortNormAstrocyteGenesZ=t(apply(ySortNormAstrocyteGenes,1,scale)); #z-score the data


meanCTRL=apply(ySortAstrocyteGenesZ[,indBlue],1,mean)
meanAD=apply(ySortAstrocyteGenesZ[,indPurple],1,mean)
logFC1=((meanAD)-(meanCTRL))

meanCTRL=apply(ySortNormAstrocyteGenesZ[,indBlue],1,mean)
meanAD=apply(ySortNormAstrocyteGenesZ[,indPurple],1,mean)
logFC2=((meanAD)-(meanCTRL))

astrocyteLogFC=data.frame(
  logFC=c(logFC1,logFC2),
  norm=rep(c("Not Normalized","Normalized"),each=length(logFC1)))


#Microglia
ySortMicrogliaGenesZ=t(apply(ySortMicrogliaGenes,1,scale)); #z-score the data
ySortNormMicrogliaGenesZ=t(apply(ySortNormMicrogliaGenes,1,scale)); #z-score the data


meanCTRL=apply(ySortMicrogliaGenesZ[,indBlue],1,mean)
meanAD=apply(ySortMicrogliaGenesZ[,indPurple],1,mean)
logFC1=((meanAD)-(meanCTRL))

meanCTRL=apply(ySortNormMicrogliaGenesZ[,indBlue],1,mean)
meanAD=apply(ySortNormMicrogliaGenesZ[,indPurple],1,mean)
logFC2=((meanAD)-(meanCTRL))

MicrogliaLogFC=data.frame(
  logFC=c(logFC1,logFC2),
  norm=rep(c("Not Normalized","Normalized"),each=length(logFC1)))







#----Normalization Plotting---- 
#Put everything in same heatmap
ySortNormZ=t(apply(ySortNorm,1,scale)); #z-score the data
ySortZNeurons=t(apply(ySortNormNeuronGenes,1,scale)); #z-score the data
ySortZAstrocytes=t(apply(ySortNormAstrocyteGenes,1,scale)); #z-score the data
ySortZMicroglia=t(apply(ySortNormMicrogliaGenes,1,scale)); #z-score the data

#ySortZNeurons[which(is.na(ySortZNeurons))]=0; #replaces NAs due to constant columns with zeros.
#No clustering


wilcox.test(ySortNormAstrocyteGenesZ[1,indPurple], ySortNormAstrocyteGenesZ[1,indBlue], alternative = "two.sided")



# ---- GSVA on Cell Type Normalized Data ----

#Choose the data to run GSVA on.
yDataforGSVA=ySortNorm;
rownames(yDataforGSVA)=geneSortClean

#----First, Wilcoxon test Astrocytes----
wilcoxTest=matrix(NA,nrow = nrow(yDataforGSVA),ncol=6)
colnames(wilcoxTest)=c("C1","SEM_C1","AD2","SEM_AD2","p","pFDR")
rownames(wilcoxTest)=rownames(yDataforGSVA)
for (i in 1:nrow(yDataforGSVA)) {
  wilcoxTest[i,1]=mean(yDataforGSVA[i,indBlue])
  wilcoxTest[i,2]=std.error(yDataforGSVA[i,indBlue])
  wilcoxTest[i,3]=mean(yDataforGSVA[i,indPurple])
  wilcoxTest[i,4]=std.error(yDataforGSVA[i,indPurple])
  wilcoxTest[i,5]=wilcox.test(yDataforGSVA[i,indPurple], yDataforGSVA[i,indBlue], alternative = "two.sided")$p.value
  
}

wilcoxTestAstrocytes=wilcoxTest[clusterNumberCTDClean==4,]
wilcoxTestAstrocytes[,6]=p.adjust(wilcoxTestAstrocytes[,5], method="fdr")
write.xlsx(wilcoxTestAstrocytes, file="astrocyteWilcoxTest_ROSMAP.xlsx", sheetName="sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)


#Choose the gene set number
geneSetNumber="1"

#Choose whether or not to compute p values
compPermute=1




#____Neurons####
geneSetsMat=read_excel(paste("../Functional Categorization/GeneSets_Neurons_FXN",geneSetNumber,".xlsx",sep=''), sheet="sheet1")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])

indclusterNumberCTDClean=which(clusterNumberCTDClean==13)


geneSetEnrich=gsva(yDataforGSVA[indclusterNumberCTDClean,], geneSets, mx.diff=TRUE, kcdf="Gaussian")

geneSetEnrich_Neurons=geneSetEnrich

#----GSVA vs Pathology----

#Braak

# indControl <-  sort(c(indBlue, indGreen))
# indMCI <- sort(c(indPlum, indIndianred))
# indAD <- sort(c(indRed, indPurple))


pathBlue <- pathScores[1,indBlue]

pathGreen <- pathScores[1,indGreen]
pathRed<- pathScores[1,indRed]
pathPurple<- pathScores[1,indPurple]
pathPlum<- pathScores[1,indPlum]
pathIndianred<- pathScores[1,indIndianred]

# indBlue.low <- colnames(pathBlue[1,which(pathBlue <4)])
# indBlue.high <- colnames(pathBlue[1,which(pathBlue >3)])
# indGreen.low <- colnames(pathGreen[1,which(pathGreen <4)])
# indGreen.high <- colnames(pathGreen[1,which(pathGreen >3)])
# indRed.low <- colnames(pathRed[1,which(pathRed <4)])
# indRed.high <- colnames(pathRed[1,which(pathRed > 3)])
# indPurple.low <- colnames(pathPurple[1,which(pathPurple <4)])
# indPurple.high <- colnames(pathPurple[1,which(pathPurple >3)])
# indPlum.low <- colnames(pathPlum[1,which(pathPlum <4)])
# indPlum.high <- colnames(pathPlum[1,which(pathPlum >3)])
# indIndianred.low <- colnames(pathIndianred[1,which(pathIndianred <4)])
# indIndianred.high <- colnames(pathIndianred[1,which(pathIndianred >3)])

pathAnalysis <- cbind(pathBlue, pathGreen, pathPurple)
indAnalysis <- colnames(pathAnalysis)
                        
indC1 <- colnames(pathBlue)                   
indC2 <- colnames(pathGreen)
indAD2 <- colnames(pathPurple)
#based on group differences in each braak score

# indPurple.0
# subgroup = indBlue.low
# folder = "AD2"
# hilo ="low"

for (i in 1:dim(geneSetEnrich)[1]){
  
  #Format data for regression
  GSVAScores=as.numeric(geneSetEnrich[13,indAnalysis])
  pathValues=as.numeric(pathScores[1,indAnalysis])
  funct.cat.name <- rownames(geneSetEnrich)
  
  pathModel=lm(GSVAScores ~ pathValues)
  
  #include number of scores in each group 
  dataPath=data.frame(pathValues,GSVAScores)
  rownames(dataPath) <- colnames(geneSetEnrich[,indAnalysis])
  dataPath$condition <- "NA"
  dataPath[indC1,]$condition <- "C1" 
  dataPath[indC2,]$condition <- "C2" 
  dataPath[indAD2,]$condition <- "AD2" 
  rowOrder <- c("C1", "C2", "AD2")
  
  dataPath <- dataPath %>%
    mutate(condition = factor(condition, levels = rowOrder)) %>%
    arrange(condition)
  
  group <- dataPath$group
  
  statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars=c("pathValues","condition"))
  
   pngPath=paste("./Braak_highlow_neuron/",folder,"/GSVA_v_Braak_",funct.cat.name[i],"_",hilo,".png",sep='')
   png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
   print({
    
    e=ggplot(data=dataPath, aes(x = as.factor(condition), y = GSVAScores, ))
    e + 
      facet_wrap(~pathValues)+
      ylim(-0.75, 0.75)+ 
      geom_violin()+
      geom_errorbar(data=statsOut, aes(group=pathValues,x= condition, ymin=GSVAScores-se, ymax=GSVAScores+se),
                    width=0.4, colour="black", alpha=1, size=1)+
      geom_crossbar(data=statsOut, aes(group=pathValues,x=condition, ymin = GSVAScores, ymax = GSVAScores),
                    size=0.5,col="red", width = 0.5)+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","AD2")))+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","C1")))+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 16),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 16),
            plot.title = element_text(hjust = 0.5))+
      xlab("")+
      ylab("Enrichment Score (a.u.)")
    
    
    
      # ggtitle(row.names(geneSetEnrich)[i])+
      # geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
      annotate("text", x = -Inf, y= -Inf, label = (paste("p = ",round(summary(pathModel)$coefficients[2,4],4))), size = 6, fontface ="bold",  hjust = -0.1, vjust= -0.5)
    

  })
  dev.off()
}



######

##### alternative version #####
# e=ggplot(data=dataPath, aes(fill = group, x = as.factor(pathValues), y = GSVAScores, ))
# e +
#   ylim(-0.55, 0.55)+
#   geom_violin(width=0.6)+
#   geom_errorbar(data=statsOut, aes(group=group,x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
#                 width=0.4, colour="black", alpha=1, size=1, position = position_dodge(width = 0.6))+
#   geom_crossbar(data=statsOut, aes(group=group,x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
#                 size=0.5,col="red", width = 0.5, position = position_dodge(width = 0.6))+
#   theme(panel.background = element_rect(),
#         text = element_text(size=20),
#         axis.text.x = element_text(face = "bold", color = "black",
#                                    size = 16),
#         axis.text.y = element_text(face = "bold", color = "black",
#                                    size = 16),
#         plot.title = element_text(hjust = 0.5))+
#   stat_compare_means(method = "wilcox.test", comparisons = list(c("C2","C1"), c("C2","AD2")))+
#   xlab("")+
#   ylab("Enrichment Score (a.u.)")+
#   ggtitle(row.names(geneSetEnrich)[i])+
#   geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
#######





####
pathBlue <- pathScores[2,indBlue]
pathGreen <- pathScores[2,indGreen]
pathRed<- pathScores[2,indRed]
pathPurple<- pathScores[2,indPurple]
pathPlum<- pathScores[2,indPlum]
pathIndianred<- pathScores[2,indIndianred]

# indBlue.low <- colnames(pathBlue[1,which(pathBlue <2)])
# indBlue.high <- colnames(pathBlue[1,which(pathBlue >1)])
# indGreen.low <- colnames(pathGreen[1,which(pathGreen <2)])
# indGreen.high <- colnames(pathGreen[1,which(pathGreen >1)])
# indRed.low <- colnames(pathRed[1,which(pathRed <2)])
# indRed.high <- colnames(pathRed[1,which(pathRed > 1)])
# indPurple.low <- colnames(pathPurple[1,which(pathPurple <2)])
# indPurple.high <- colnames(pathPurple[1,which(pathPurple >1)])
# indPlum.low <- colnames(pathPlum[1,which(pathPlum <2)])
# indPlum.high <- colnames(pathPlum[1,which(pathPlum >1)])
# indIndianred.low <- colnames(pathIndianred[1,which(pathIndianred <2)])
# indIndianred.high <- colnames(pathIndianred[1,which(pathIndianred >1)])

pathAnalysis <- cbind(pathBlue, pathGreen, pathPurple)
indAnalysis <- colnames(pathAnalysis)

indC1 <- colnames(pathBlue)                   
indC2 <- colnames(pathGreen)
indAD2 <- colnames(pathPurple)

# subgroup = indPurple.high
# folder = "AD2"
# hilo ="high"


#CERAD
for (i in 1:dim(geneSetEnrich)[1]){
  #Format data for regression
  #Format data for regression
  GSVAScores=as.numeric(geneSetEnrich[1,indAnalysis])
  pathValues=as.numeric(pathScores[2,indAnalysis])
  funct.cat.name <- rownames(geneSetEnrich)
  
  pathModel=lm(GSVAScores ~ pathValues)
  
  #include number of scores in each group 
  dataPath=data.frame(pathValues,GSVAScores)
  rownames(dataPath) <- colnames(geneSetEnrich[,indAnalysis])
  dataPath$condition <- "NA"
  dataPath[indC1,]$condition <- "C1" 
  dataPath[indC2,]$condition <- "C2" 
  dataPath[indAD2,]$condition <- "AD2" 
  rowOrder <- c("C1", "C2", "AD2")
  
  dataPath <- dataPath %>%
    mutate(condition = factor(condition, levels = rowOrder)) %>%
    arrange(condition)
  
  group <- dataPath$group
  
  statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars=c("pathValues","condition"))
  
  pngPath=paste("./Braak_highlow_neuron/",folder,"/GSVA_v_Braak_",funct.cat.name[i],"_",hilo,".png",sep='')
  png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
  print({
    
    e=ggplot(data=dataPath, aes(x = as.factor(condition), y = GSVAScores, ))
    e + 
      facet_wrap(~pathValues)+
      ylim(-0.7, 0.7)+ 
      geom_violin()+
      geom_errorbar(data=statsOut, aes(group=pathValues,x= condition, ymin=GSVAScores-se, ymax=GSVAScores+se),
                    width=0.4, colour="black", alpha=1, size=1)+
      geom_crossbar(data=statsOut, aes(group=pathValues,x=condition, ymin = GSVAScores, ymax = GSVAScores),
                    size=0.5,col="red", width = 0.5)+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","AD2")))+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","C1")))+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 16),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 16),
            plot.title = element_text(hjust = 0.5))+
      xlab("")+
      ylab("Enrichment Score (a.u.)")+
      ggtitle(row.names(geneSetEnrich)[i])+
      geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
  
  })
  dev.off()
  
}



#____Astrocytes####
geneSetsMat=read_excel(paste("../Functional Categorization/GeneSets_Astrocytes_FXN", geneSetNumber, ".xlsx",sep=''), sheet="sheet1")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])

geneSetEnrich=gsva(yDataforGSVA[clusterNumberCTDClean==4,], geneSets, mx.diff=TRUE, kcdf="Gaussian")

geneSetEnrich_Astrocytes=geneSetEnrich

pathBlue <- pathScores[1,indBlue]
pathGreen <- pathScores[1,indGreen]
pathRed<- pathScores[1,indRed]
pathPurple<- pathScores[1,indPurple]
pathPlum<- pathScores[1,indPlum]
pathIndianred<- pathScores[1,indIndianred]

# indBlue.low <- colnames(pathBlue[1,which(pathBlue <4)])
# indBlue.high <- colnames(pathBlue[1,which(pathBlue >3)])
# indGreen.low <- colnames(pathGreen[1,which(pathGreen <4)])
# indGreen.high <- colnames(pathGreen[1,which(pathGreen >3)])
# indRed.low <- colnames(pathRed[1,which(pathRed <4)])
# indRed.high <- colnames(pathRed[1,which(pathRed > 3)])
# indPurple.low <- colnames(pathPurple[1,which(pathPurple <4)])
# indPurple.high <- colnames(pathPurple[1,which(pathPurple >3)])
# indPlum.low <- colnames(pathPlum[1,which(pathPlum <4)])
# indPlum.high <- colnames(pathPlum[1,which(pathPlum >3)])
# indIndianred.low <- colnames(pathIndianred[1,which(pathIndianred <4)])
# indIndianred.high <- colnames(pathIndianred[1,which(pathIndianred >3)])

pathAnalysis <- cbind(pathBlue, pathGreen, pathPurple)
indAnalysis <- colnames(pathAnalysis)

indC1 <- colnames(pathBlue)
indC2 <- colnames(pathGreen)
indAD2 <- colnames(pathPurple)

#----GSVA vs Pathology----

# subgroup = indPurple.high
# folder = "AD2"
# hilo ="high"

#Braak
pathModel_pVal=matrix(NA,nrow=dim(geneSetEnrich)[1], ncol = 1)

for (i in 1:dim(geneSetEnrich)[1]){
  
  #Format data for regression
  GSVAScores=as.numeric(geneSetEnrich[14,indAnalysis])
  pathValues=as.numeric(pathScores[1,indAnalysis])
  funct.cat.name <- rownames(geneSetEnrich)
  
  pathModel=lm(GSVAScores ~ pathValues)
  
  #include number of scores in each group 
  dataPath=data.frame(pathValues,GSVAScores)
  rownames(dataPath) <- colnames(geneSetEnrich[,indAnalysis])
  dataPath$condition <- "NA"
  dataPath[indC1,]$condition <- "C1" 
  dataPath[indC2,]$condition <- "C2" 
  dataPath[indAD2,]$condition <- "AD2" 
  rowOrder <- c("C1", "C2", "AD2")
  
  dataPath <- dataPath %>%
    mutate(condition = factor(condition, levels = rowOrder)) %>%
    arrange(condition)
  
  group <- dataPath$group
  
  statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars=c("pathValues","condition"))
  

  pngPath=paste("./Braak_highlow_astrocyte/",folder,"/GSVA_v_Braak_",funct.cat.name[i],"_",hilo,".png",sep='')
  png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
  print({
    
    e=ggplot(data=dataPath, aes(x = as.factor(condition), y = GSVAScores, ))
    e + 
      facet_wrap(~pathValues)+
      ylim(-0.8, 0.8)+ 
      geom_violin()+
      geom_errorbar(data=statsOut, aes(group=pathValues,x= condition, ymin=GSVAScores-se, ymax=GSVAScores+se),
                    width=0.4, colour="black", alpha=1, size=1)+
      geom_crossbar(data=statsOut, aes(group=pathValues,x=condition, ymin = GSVAScores, ymax = GSVAScores),
                    size=0.5,col="red", width = 0.5)+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","AD2")))+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","C1")))+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 16),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 16),
            plot.title = element_text(hjust = 0.5))+
      xlab("")+
      ylab("Enrichment Score (a.u.)")+
      
      ggtitle(row.names(geneSetEnrich)[i])+
      geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
      # 
     })
  dev.off()
}

# pathModelOut=data.frame(geneSet=rownames(geneSetEnrich),pValSlope=pathModel_pVal)
# write.xlsx(pathModelOut, file="BraakvsAstGeneSets_MCI2.xlsx", sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)


pathBlue <- pathScores[2,indBlue]
pathGreen <- pathScores[2,indGreen]
pathRed<- pathScores[2,indRed]
pathPurple<- pathScores[2,indPurple]
pathPlum<- pathScores[2,indPlum]
pathIndianred<- pathScores[2,indIndianred]

# indBlue.low <- colnames(pathBlue[1,which(pathBlue <2)])
# indBlue.high <- colnames(pathBlue[1,which(pathBlue >1)])
# indGreen.low <- colnames(pathGreen[1,which(pathGreen <2)])
# indGreen.high <- colnames(pathGreen[1,which(pathGreen >1)])
# indRed.low <- colnames(pathRed[1,which(pathRed <2)])
# indRed.high <- colnames(pathRed[1,which(pathRed > 1)])
# indPurple.low <- colnames(pathPurple[1,which(pathPurple <2)])
# indPurple.high <- colnames(pathPurple[1,which(pathPurple >1)])
# indPlum.low <- colnames(pathPlum[1,which(pathPlum <2)])
# indPlum.high <- colnames(pathPlum[1,which(pathPlum >1)])
# indIndianred.low <- colnames(pathIndianred[1,which(pathIndianred <2)])
# indIndianred.high <- colnames(pathIndianred[1,which(pathIndianred >1)])

pathAnalysis <- cbind(pathBlue, pathGreen, pathPurple)
indAnalysis <- colnames(pathAnalysis)

indC1 <- colnames(pathBlue)
indC2 <- colnames(pathGreen)
indAD2 <- colnames(pathPurple)
# 
# subgroup = indPurple.high
# folder = "AD2"
# hilo ="high"


#CERAD
pathModel_pVal=matrix(NA,nrow=dim(geneSetEnrich)[1], ncol = 1)

for (i in 1:dim(geneSetEnrich)[1]){
  
  #Format data for regression
  GSVAScores=as.numeric(geneSetEnrich[14,indAnalysis])
  pathValues=as.numeric(pathScores[2,indAnalysis])
  funct.cat.name <- rownames(geneSetEnrich)
  
  pathModel=lm(GSVAScores ~ pathValues)
  
  #include number of scores in each group 
  dataPath=data.frame(pathValues,GSVAScores)
  rownames(dataPath) <- colnames(geneSetEnrich[,indAnalysis])
  dataPath$condition <- "NA"
  dataPath[indC1,]$condition <- "C1" 
  dataPath[indC2,]$condition <- "C2" 
  dataPath[indAD2,]$condition <- "AD2" 
  rowOrder <- c("C1", "C2", "AD2")
  
  dataPath <- dataPath %>%
    mutate(condition = factor(condition, levels = rowOrder)) %>%
    arrange(condition)
  
  group <- dataPath$group
  
  statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars=c("pathValues","condition"))
  
  
  pngPath=paste("./CERAD_highlow_astrocyte/",folder,"/GSVA_v_CERAD_",funct.cat.name[i],"_",hilo,".png",sep='')
  png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
  print({
    
    e=ggplot(data=dataPath, aes(x = as.factor(condition), y = GSVAScores, ))
    e + 
      facet_wrap(~pathValues)+
      ylim(-0.7, 0.7)+ 
      geom_violin()+
      geom_errorbar(data=statsOut, aes(group=pathValues,x= condition, ymin=GSVAScores-se, ymax=GSVAScores+se),
                    width=0.4, colour="black", alpha=1, size=1)+
      geom_crossbar(data=statsOut, aes(group=pathValues,x=condition, ymin = GSVAScores, ymax = GSVAScores),
                    size=0.5,col="red", width = 0.5)+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","AD2")))+
      stat_compare_means(method="wilcox.test", comparisons = list(c("C2","C1")))+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 16),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 16),
            plot.title = element_text(hjust = 0.5))+
      xlab("")+
      ylab("Enrichment Score (a.u.)")+
      
      ggtitle(row.names(geneSetEnrich)[i])+
      geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
    # 
  })
  dev.off()
}
  

# pathModelOut=data.frame(geneSet=rownames(geneSetEnrich),pValSlope=pathModel_pVal)
# write.xlsx(pathModelOut, file="CERADvsAstGeneSets_AD2.xlsx", sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)
