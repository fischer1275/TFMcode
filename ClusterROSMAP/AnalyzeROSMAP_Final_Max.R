# ---- Preliminaries ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

pngPath=paste("Normalization__neurons_hist.png",sep='')
png(pngPath, width=4,height=4,units="in",res=600, pointsize = 20)
ggdensity(neuoronLogFC, x="logFC", y = "..density..",
          color="norm", fill = "norm", palette = "jco",
          size=1, )+labs(title="Neuron Genes")+xlab("Mean(AD2-C1)")+ylab("Density")+xlim(-3,3)+
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1.8))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        text = element_text(size=12),
        axis.text.x = element_text(face = "plain", color = "black", 
                                   size = 24),
        axis.title.x = element_text(face = "plain", color = "black", 
                                    size = 30),
        axis.text.y = element_text(face = "plain", color = "black", 
                                   size = 24),
        axis.title.y = element_text(face = "plain", color = "black", 
                                    size = 30))
dev.off()

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

png("Normalization_astrocytes_hist.png", width=4,height=4,units="in",res=600, pointsize = 20)

ggdensity(astrocyteLogFC, x="logFC", y = "..density..",
          color="norm", fill = "norm", palette = "jco",
          size=1, )+labs(title="Astrocyte Genes")+xlab("Mean (A2-C1)")+ylab("Density")+xlim(-3,3)+ylim(0,1.3)+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        text = element_text(size=12),
        axis.text.x = element_text(face = "plain", color = "black", 
                                   size = 24),
        axis.title.x = element_text(face = "plain", color = "black", 
                                    size = 30),
        axis.text.y = element_text(face = "plain", color = "black", 
                                   size = 24),
        axis.title.y = element_text(face = "plain", color = "black", 
                                    size = 30))
dev.off()

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

png("Normalization_Microglias_hist.png", width=4,height=4,units="in",res=600, pointsize = 20)

ggdensity(MicrogliaLogFC, x="logFC", y = "..density..",
          color="norm", fill = "norm", palette = "jco",
          size=1, )+labs(title="Microglia Genes")+xlab("Mean (A2-C1)")+ylab("Density")+xlim(-3,3)+ylim(0,1.3)+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'), 
        text = element_text(size=12),
        axis.text.x = element_text(face = "plain", color = "black", 
                                   size = 24),
        axis.title.x = element_text(face = "plain", color = "black", 
                                    size = 30),
        axis.text.y = element_text(face = "plain", color = "black", 
                                   size = 24),
        axis.title.y = element_text(face = "plain", color = "black", 
                                    size = 30))
dev.off()




#@@@@@ extract cell type specific dendrograms  @@@@@·

# genes.cluster.pair <- data.frame(geneLabelsClustCTD, clusterNumberCTD) 
# 
# write.xlsx(genes.cluster.pair, "./genenames_cluster_pair.xlsx")
# 
# #gene names NOT matched with with database genes 
# clust4.genes <- genes.cluster.pair[genes.cluster.pair[,2]==4,1] 
# clust13.genes <- genes.cluster.pair[genes.cluster.pair[,2]==13,1] 
# clust17.genes <- genes.cluster.pair[genes.cluster.pair[,2]==17,1]  
# 
# hrCTD.in <- yMatCTDZ
# rownames(hrCTD.in) <- geneLabelsCTD
# 
# hrCTD.2 = hclust(dist((hrCTD.in),method = "euclidean"), method="average")  
# 
# sub_dent.4 <- find_dendrogram(as.dendrogram(hrCTD.2), clust4.genes)
# 
# sub_dent.13 <- find_dendrogram(as.dendrogram(hrCTD.2), clust13.genes)
# 
# sub_dent.17 <- find_dendrogram(as.dendrogram(hrCTD.2), clust17.genes)
# 
# 
# png("Sub dendrogram cluster 4.png")
# plot(sub_dent.4, horiz=T, labels= F)
# dev.off()
# 
# png("Sub dendrogram cluster 13.png")
# plot(sub_dent.13, horiz=T, labels= F)
# dev.off()
# 
# png("Sub dendrogram cluster 17.png")
# plot(sub_dent.17, horiz=T, labels= F)
# dev.off()
# 
# 
# test1 <- ySortNorm
# rownames(test1) <- geneSortClean
# 
# plot.norm.neu <- test1[clusterNumberCTDClean==13,]
# plot.norm.ast  <- test1[clusterNumberCTDClean==4,]
# plot.norm.mic  <- test1[clusterNumberCTDClean==17,]
# 
# #these are now the same values as ySortZNeuron/Astro/Micro but with the gene names as row names 
# plot.norm.neu.z <- t(apply(plot.norm.neu,1,scale))
# plot.norm.neu.z <- plot.norm.neu.z[-c(1:5),]
# plot.norm.neu.z <- as.data.frame(plot.norm.neu.z)
# 
# plot.norm.ast.z  <- t(apply(plot.norm.ast,1,scale))
# plot.norm.mic.z <- t(apply(plot.norm.mic,1,scale))
# 
# 
# plot.data <- ySortZ
# rownames(plot.data) <- geneLabelsClustCTD
# 
# plot.data.neu <- plot.data[clusterNumberCTD==13,]
# plot.data.neu <- as.data.frame(plot.data.neu)
# 
# plot.data.ast <- plot.data[clusterNumberCTD==4,]
# plot.data.mic <- plot.data[clusterNumberCTD==17,]
# 
# plot.data.neu[is.na(plot.data.neu)] <- 0 # set NA values to 0 
# 
# library(tibble)
# 
# #updates the Zhang values in the present genes with the VSN normalised and z-scored data from the dataset (values are now equeal to ySortZNeurons/Astrocytes/mciroglia)
# t3 <- plot.data.neu %>%
#   rownames_to_column() %>%
#   rows_update(plot.norm.neu.z %>% rownames_to_column(), by ="rowname") %>%
#   column_to_rownames()
# 
# #plot the heatmap
# heatmap3(t3[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
#          col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
#          Rowv= sub_dent.13, Colv=NA,  scale="none",
#          labRow = NA, labCol = NA, main="Normalized Neuron Genes")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@2



#----Normalization Plotting---- 
#Put everything in same heatmap
ySortNormZ=t(apply(ySortNorm,1,scale)); #z-score the data
ySortZNeurons=t(apply(ySortNormNeuronGenes,1,scale)); #z-score the data
ySortZAstrocytes=t(apply(ySortNormAstrocyteGenes,1,scale)); #z-score the data
ySortZMicroglia=t(apply(ySortNormMicrogliaGenes,1,scale)); #z-score the data
  
#ySortZNeurons[which(is.na(ySortZNeurons))]=0; #replaces NAs due to constant columns with zeros.
#No clustering


#Column clustering whole dataset
hcSort= hclust(dist(t(ySortNormZ),method = "euclidean"), method="complete")
heatmap3(ySortNormZ, ColSideColors =c(colorBar), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSort),  scale="row",
         labRow = NA, labCol = NA) 

#hrSort = hclust(as.dist(1-cor(t(ySortZNeurons), method="pearson", use="na.or.complete")), method="complete") 


#Normalized Neuron and Astrocyte Genes (WITH subgroup filtering)
heatmap3(ySortZNeurons[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="none",
         labRow = NA, labCol = NA, main="Normalized Neuron Genes")

heatmap3(ySortZAstrocytes[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="none",
         labRow = NA, labCol = NA, main="Normalized Astrocyte Genes") 

heatmap3(ySortZMicroglia[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="none",
         labRow = NA, labCol = NA, main="Normalized Microglia Genes") 


#NOT Normalized Neuron and Astrocyte Genes for comparison
heatmap3(ySortNeuronGenes[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA, main=" Neuron Genes")

heatmap3(ySortAstrocyteGenes[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA, main="Astrocyte Genes") 

heatmap3(ySortMicrogliaGenes[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA, main="Microglia Genes") 


wilcox.test(ySortNormAstrocyteGenesZ[1,indPurple], ySortNormAstrocyteGenesZ[1,indBlue], alternative = "two.sided")



# ---- Unsupervised self-clustering of Cell Type specific gene groups --- #

###NEURONS

clust.neurons_avg <- hclust(dist((ySortZNeurons), method = "euclidean"), method ="average")
clust.neurons_ward <- hclust(dist((ySortZNeurons), method = "euclidean"), method ="ward.D2")
cort <- cor(t(ySortZNeurons), method="spearman")
# clust.neurons.2 <- hclust(as.dist(1-cort),method = "complete")
clust.neurons <- eclust(ySortZNeurons, "hclust", k=NULL, hc_metric = "spearman", hc_method = "complete", graph = FALSE )

#plot(as.dendrogram(clust.neurons))

# gap statistics for internal cluster number validation
png("Gap statistics neurons.png", width=4,height=4,units="in",res=600, pointsize = 20)

fviz_nbclust(ySortZNeurons, hcut, method = "gap_stat", diss = as.dist(1-cort), k.max = 5)

dev.off()

#cut dendrogram into indicated clusters
clust.neurons.cut <- cutree(clust.neurons_ward, k=5)

colorBrewer.clust.neu =colorRampPalette(brewer.pal(8,"Spectral"))
mycolclust.neu =  sample(colorBrewer.clust.neu(length(unique(clust.neurons.cut))))
mycolclust.neu = mycolclust.neu[as.vector(clust.neurons.cut)]

png("Norm Neuron Genes self clustering (ward) k5.png", width=4,height=4,units="in",res=600, pointsize = 20)

heatmap3(ySortZNeurons[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv= as.dendrogram(clust.neurons_ward), Colv=NA,  scale="none",
         RowSideColors = mycolclust.neu, RowSideLabs=NA,
         labRow = NA, labCol = NA, main="Normalized Neuron Genes (celltype individual clustering)")

dev.off()


### ASTROCYTES ###

clust.astrocytes_avg <- hclust(dist((ySortZAstrocytes), method = "euclidean"), method ="average")
clust.astrocytes_ward <- hclust(dist((ySortZAstrocytes), method = "euclidean"), method ="ward.D2")
cort.ast <- cor(t(ySortZAstrocytes), method = "spearman")
# clust.astrocytes <- hclust(as.dist(1-cort),method = "complete")
clust.astrocytes <- eclust(ySortZAstrocytes, "hclust", k=NULL, hc_metric = "spearman", hc_method = "complete", graph = FALSE )

#gap statistics for cluster number validation
png("Gap statistics Astrocytes.png", width=4,height=4,units="in",res=600, pointsize = 20)

fviz_nbclust(ySortZAstrocytes, hcut, method = "gap_stat", diss = as.dist(1-cort.ast), k.max = 5)

dev.off()

#plot(as.dendrogram(clust.astrocytes))

clust.astrocytes.cut <- cutree(clust.astrocytes_ward, k = 5)

colorBrewer.clust.ast =colorRampPalette(brewer.pal(8,"Spectral"))
mycolclust.ast =  sample(colorBrewer.clust.ast(length(unique(clust.astrocytes.cut))))
mycolclust.ast = mycolclust.ast[as.vector(clust.astrocytes.cut)]

png("Norm Astrocyte Genes self clustering (ward) k5.png", width=4,height=4,units="in",res=600, pointsize = 20)

heatmap3(ySortZAstrocytes[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv= as.dendrogram(clust.astrocytes_ward), Colv=NA,  scale="none",
         RowSideColors = mycolclust.ast, RowSideLabs=NA,
         labRow = NA, labCol = NA, main="Normalized Astrocyte Genes (celltype individual clustering)")

dev.off()

### MICROGLIA ###

# there are 3 genes that are 0 but get converted to NA, need to convert tehm back to 0 
ySortZMicroglia[is.na(ySortZMicroglia)] <- 0

ySortZMicroglia_ann <- ySortZMicroglia
#rownames(ySortZMicroglia_ann) <- microGeneNames

# common <- read.xlsx("../Functional Categorization/CommonMicrogliaGenes.xlsx", sheetIndex = 1, header = FALSE)
# common_list <- as.list(common)
# common_keep <- common_list$X1
# 
# yMicrogliaClean <- subset(ySortZMicroglia_ann, rownames(ySortZMicroglia_ann) %in% common_keep)

clust.microglia_avg <- hclust(dist((ySortZMicroglia), method = "euclidean"), method ="average")
clust.microglia_ward <- hclust(dist((ySortZMicroglia), method = "euclidean"), method ="ward.D2")
cort.micro <- cor(t(ySortZMicroglia), method="spearman")
# clust.microglia <- hclust(as.dist(1-cort),method = "complete")

clust.microglia <- eclust(ySortZMicroglia, "hclust", k=NULL, hc_metric = "spearman", hc_method = "complete", graph = FALSE )

#gap statistics for cluster number validation
png("Gap statistics microglia.png", width=4,height=4,units="in",res=600, pointsize = 20)

fviz_nbclust(ySortZMicroglia, hcut, method = "gap_stat", diss = as.dist(1-cort.micro), k.max = 5)

dev.off()

#plot(as.dendrogram(clust.microglia))

clust.microglia.cut <- cutree(clust.microglia, k=5) 

#table(clust.microglia.cut)

#cluster numbers as shown in heatmap
# A - clust 5 (164)
# B - clust 3 (109)
# C - clust 4 (280)
# D - clust 2 (288)
# E - clust 1 (325)

colorBrewer.clust.mic =colorRampPalette(brewer.pal(8,"Spectral"))
mycolclust.mic =  sample(colorBrewer.clust.mic(length(unique(clust.microglia.cut))))
mycolclust.mic = mycolclust.mic[as.vector(clust.microglia.cut)]

png("Norm Microglia Genes self clustering (corr) k5.png", width=4,height=4,units="in",res=600, pointsize = 20)

heatmap3(ySortZMicroglia[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv= as.dendrogram(clust.microglia), Colv=NA,  scale="none",
         RowSideColors = mycolclust.mic, RowSideLabs=NA,
         labRow = NA, labCol = NA, main="Normalized Microglia Genes (celltype individual clustering)")

dev.off()


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
    
    indBlue.low <- colnames(pathBlue[1,which(pathBlue <3)])
    indBlue.high <- colnames(pathBlue[1,which(pathBlue >2)])
    indGreen.low <- colnames(pathGreen[1,which(pathGreen <3)])
    indGreen.high <- colnames(pathGreen[1,which(pathGreen >2)])
    indRed.low <- colnames(pathRed[1,which(pathRed <3)])
    indRed.high <- colnames(pathRed[1,which(pathRed > 2)])
    indPurple.low <- colnames(pathPurple[1,which(pathPurple <3)])
    indPurple.high <- colnames(pathPurple[1,which(pathPurple >2)])
    indPlum.low <- colnames(pathPlum[1,which(pathPlum <3)])
    indPlum.high <- colnames(pathPlum[1,which(pathPlum >2)])
    indIndianred.low <- colnames(pathIndianred[1,which(pathIndianred <3)])
    indIndianred.high <- colnames(pathIndianred[1,which(pathIndianred >2)])
    
    
    for (i in 1:dim(geneSetEnrich)[1]){
      
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,indBlue.low])
      pathValues=as.numeric(pathScores[1,indBlue.low])
    
      funct.cat.name <- rownames(geneSetEnrich)
      
      pathModel=lm(GSVAScores ~ pathValues)
      
      #include number of scores in each group 
      
      
      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
      
      pngPath=paste("./Braak_highlow_neuron/C1/GSVA_v_Braak_",funct.cat.name[i],"_low.png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({
        
        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1)+
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "bold", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "bold", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("Braak Score")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])+
          geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])+
          annotate("text", x = -Inf, y= -Inf, label = (paste("p = ",round(summary(pathModel)$coefficients[2,4],4))), size = 6, fontface ="bold",  hjust = -0.1, vjust= -0.5)
          
        #paste("p =",summary(pathModel)$coefficients[2,4])
      })
      dev.off()
      
    }
    
    
    #CERAD
    for (i in 1:dim(geneSetEnrich)[1]){
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,indBlue.low])
      pathValues=as.numeric(pathScores[2,indBlue.low])
      pathModel=lm(GSVAScores ~ pathValues)
      
      funct.cat.name <- rownames(geneSetEnrich)
      
      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
      
      pngPath=paste("./CERAD_highlow_neuron/C1/GSVA_vs_CERAD_",funct.cat.name[i],"_low.png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({
        
        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1) +
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "plain", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("CERAD")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])+
          geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
      })
      dev.off()
      
    }
    
    #Plot CERAD vs Braak
    x=as.matrix(pathScores[1,])
    y=as.matrix(pathScores[2,])

      
    plot(x,y, xlab="Braak", 
         ylab="CERAD")
    
    
    #APOE
    for (i in 1:dim(geneSetEnrich)[1]){
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,])
      pathValues=as.numeric(pathScores[3,])

      GSVAScores=GSVAScores[-which(is.na(pathValues))]
      pathValues=pathValues[-which(is.na(pathValues))]

      pathModel=lm(GSVAScores ~ pathValues)

      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")

      pngPath=paste("./FXN1_norm_Neurons/GSVA_vs_APOE_",i,".png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({

        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1) +
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "plain", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("APOE")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])
      })
      dev.off()

    }

    #Sex
    for (i in 1:dim(geneSetEnrich)[1]){
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,])
      pathValues=as.numeric(pathScores[4,])
      pathModel=lm(GSVAScores ~ pathValues)

      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")

      pngPath=paste("./FXN1_norm_Neurons/GSVA_vs_Sex_",i,".png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({

        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1) +
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "plain", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("Sex")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])+
          scale_x_discrete(labels=c("0" = "Female", "1" = "Male"))+
          geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
      })
      dev.off()
      
    }
    #----GSVA heatmap----

    
    geneSetEnrichSort=geneSetEnrich_Neurons[,indSort] #sort samples by phenotype
    geneSetEnrichSortZ=t(apply(geneSetEnrich_Neurons,1,scale)); #z-score the data
    
    hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons/GSVA_cluster.png",sep='')
    png(pngPath, width=8,height=4,units="in",res=600, pointsize = 16)

    heatmap3(geneSetEnrichSortZ, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
             col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
             Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
             labRow = colnames(geneSetsMat), labCol = NA, cexRow=1.2) 
    title("GSVA for neuron gene sets", adj = 0.5, line = 2)
    
    dev.off()
    
    
    corGSVA=cor(t(geneSetEnrichSortZ))
    hCor=hclust(as.dist((1-corGSVA)/2))
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons/GSVA_corr.png",sep='')
    png(pngPath, width=6,height=8,units="in",res=600, pointsize = 14)
    #par(mar=c(5,6,4,4)+.1)
    
    heatmap3(corGSVA,
             col=barColorsCor, breaks=breakBarColorsCor,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
             Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",
             labRow = colnames(geneSetsMat), labCol = colnames(geneSetsMat),margins=c(12,13)) 
      title("GSVA for neuron gene sets", adj = 0.5, line = 1)
    dev.off()
    
    blah=yDataforGSVA[indclusterNumberCTDClean,]
    blah_C2AD1 <- 
    
    #Permutation p-value
    if (compPermute==1){
      R=1000
      meanOut6=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut5=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut4=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut3=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut2=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut1=matrix(NA,nrow=length(geneSets), ncol = R)
      
      for (i in 1:R)
      {
        blah2=blah
        rownames(blah2)=rownames(blah)[sample(1:dim(blah)[1])]
        gsvaBoot=gsva(blah2,geneSets, mx.diff=TRUE, kcdf="Gaussian")
        meanOut6[,i]=rowMeans(gsvaBoot[,indIndianred])-rowMeans(gsvaBoot[,indPlum]) # MCI2 - MCI1
        meanOut5[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPurple]) # AD2 - C2
        meanOut4[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indRed]) # AD1 - C2 
        meanOut3[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indIndianred]) # MCI2 - C2
        meanOut2[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPlum]) #MCI1 - C2
        meanOut1[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indBlue]) # C1 - C2
        
        
        print(i)
      }
      meanTrue6=rowMeans(geneSetEnrich[,indIndianred])-rowMeans(geneSetEnrich[,indPlum])
      meanTrue5=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPurple])
      meanTrue4=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indRed])
      meanTrue3=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indIndianred])
      meanTrue2=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPlum])
      meanTrue1=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indBlue])
      
      
      pBoot6=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot5=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot4=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot3=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot2=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot1=matrix(NA,nrow=length(geneSets), ncol = 1)
      
      for (j in 1:length(geneSets))
        
      {
        pBoot6[j] = mean(abs(meanOut6[j,]) > abs(meanTrue6[j]))
        pBoot5[j] = mean(abs(meanOut5[j,]) > abs(meanTrue5[j]))
        pBoot4[j] = mean(abs(meanOut4[j,]) > abs(meanTrue4[j]))
        pBoot3[j] = mean(abs(meanOut3[j,]) > abs(meanTrue3[j]))
        pBoot2[j] = mean(abs(meanOut2[j,]) > abs(meanTrue2[j]))
        pBoot1[j] = mean(abs(meanOut1[j,]) > abs(meanTrue1[j]))
        
      }
      pBootFDR6=p.adjust(pBoot6,method = "fdr",)
      pBootFDR5=p.adjust(pBoot5,method = "fdr",)
      pBootFDR4=p.adjust(pBoot4,method = "fdr",)
      pBootFDR3=p.adjust(pBoot3,method = "fdr",)
      pBootFDR2=p.adjust(pBoot2,method = "fdr",)
      pBootFDR1=p.adjust(pBoot1,method = "fdr",)
      
      
      pBoot=data.frame(pBoot1, pBoot2, pBoot3, pBoot4, pBoot5, pBoot6)
      pBootFDR=data.frame(pBootFDR1, pBootFDR2, pBootFDR3, pBootFDR4, pBootFDR5, pBootFDR6)
      
      rownames(pBoot)=rownames(geneSetEnrich)
      rownames(pBootFDR)=rownames(geneSetEnrich)
      #End permutation
      
      pBootNeu=pBoot
      pBootFDRNeu=pBootFDR
      
    }
    
    if(geneSetNumber=="1"){
      save_pBootNeu1=data.frame(pBootNeu,pBootFDRNeu)
    }else{
      save_pBootNeu2=data.frame(pBootNeu,pBootFDRNeu)
    }
    
    
    #plotting neuron genes
    heatmap3(yDataforGSVA[clusterNumberCTDClean==13,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
             col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
             Rowv=NA, Colv=NA,  scale="row",
             labRow = NA, labCol = NA, main="neuron genes") 
    
    #Compute p-values and FDR q-values for each gene set.
    # pValsTukey=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
    # pANOVA=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
    numGenes=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
    
    for (i in 1:dim(geneSetEnrich)[1])
    { #p=t.test(geneSetEnrichSort[i,allColorGroups=="Blue"], geneSetEnrichSort[i,allColorGroups=="Purple"])
      #pVals[i]=p$p.value
      
      indGetGS=which(match(names(geneSets),row.names(geneSetEnrich)[i])==1)
      numGenes[i]=length(unlist(geneSets[indGetGS]))
      
      #ANOVA
    #   dataANOVA=data.frame(geneSetEnrichSort[i,], allColorGroups)
    #   colnames(dataANOVA)=c("GSVA_ES","color")
    #   AOVObject=aov(GSVA_ES ~ color, data = dataANOVA)
    #   resultANOVA=summary(AOVObject)
    #   pANOVA[i]=resultANOVA[[1]][["Pr(>F)"]][1]
    #   
    #   ANOVATukey=TukeyHSD(AOVObject,conf.level = 0.95)
    #   pValsTukey[i]=ANOVATukey$color[which(rownames(ANOVATukey$color)=="skyblue4-orangered4"),4] #Compare blue vs purple
     }
    
    # FDRANOVA=p.adjust(pANOVA,method = "fdr")
    # qFDRTukey=p.adjust(pValsTukey,method = "fdr")
    geneSets_Stats=data.frame(rownames(geneSetEnrichSort),pBootNeu,pBootFDRNeu, numGenes);
    write.xlsx(geneSets_Stats, file=paste("geneSetStats_Neurons_FXN",geneSetNumber,"_norm_C2vALL.xlsx",sep=''), sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)
    
    
    
    #Plot bar graphs and gene set heatmaps
    
    meanBluevsOR=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) #C2 v AD1
    meanBluevsOR4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v AD2
    meanBluevsPapaya=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v MCI1
    meanBluevsPeach=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v MCI2
    meanBluevsBlue4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v C1
    meanPeachvsPapaya =matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) #MCI2 v MCI1
    
    for (indGeneSetPlot in 1:dim(geneSetEnrich)[1])
    {
      dataBar=data.frame(as.data.frame(geneSetEnrich[indGeneSetPlot,indSort]), as.data.frame(allColorGroups))
      colnames(dataBar)=c("ES","colors")
      
      #Summarize stats for each group color
      statsOut = summarySE(dataBar, measurevar="ES", groupvars="colors")
      
      meanBluevsOR[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered"] # AD1 - C2
      meanBluevsOR4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered4"] #AD2 - C2
      meanBluevsPapaya[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="papayawhip"] #C2 v MC1
      meanBluevsPeach[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="peachpuff3"] #C2 vs MC2
      meanBluevsBlue4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="skyblue4"] # C2 vs C1
      meanPeachvsPapaya[indGeneSetPlot]=statsOut$ES[statsOut$color=="peachpuff3"]-statsOut$ES[statsOut$color=="papayawhip"] # MCI2 v MCI1
      
      
    
      #Re-order the colors to match the heatmaps
      statsOut$colors = factor(statsOut$colors, levels = c("skyblue4","skyblue","papayawhip","peachpuff3","orangered","orangered4"))
      
      group.colors <- c(skyblue4 = "skyblue4", skyblue = "skyblue", papayawhip="papayawhip",peachpuff3="peachpuff3", orangered ="orangered", orangered4 = "orangered4")
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons/geneSets",indGeneSetPlot,".png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600)
      print({
        ggplot(statsOut, aes(x=colors, y=ES, fill=colors)) + 
            geom_bar(position=position_dodge(), stat="identity", colour="black") +
            geom_errorbar(aes(ymin=ES-se, ymax=ES+se),
                          width=.2, size=1.5,                    # Width of the error bars
                          position=position_dodge(.9))+
            theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                  text = element_text(size=28),
                  axis.text.x = element_text(face = "plain", color = "black", 
                                             size = 24),
                  axis.text.y = element_text(face = "plain", color = "black", 
                                             size = 24),
                  legend.position = "none")+
          xlab("")+
          ylab("Enrichment Score")+
            ggtitle(rownames(geneSetEnrich)[indGeneSetPlot])+
            scale_fill_manual(values=group.colors)+
            ylim(-0.5,0.5)+
          scale_x_discrete(labels=c("C1", "C2", "M1", "M2", "A1", "A2"))+
          geom_hline(yintercept=0)
        
      })
      dev.off()
      
      #plot heatmaps for each gene set
      
        #Get index of gene set associated with enrichment
        indGetGS=which(match(names(geneSets),row.names(geneSetEnrich)[indGeneSetPlot])==1)
        
      
      indGS_inData=match(as.vector(unlist(geneSets[indGetGS])),geneSortClean)
      indGS_inData=indGS_inData[!is.na(indGS_inData)]
      
      FXN_Annotation_Neurons=read_excel("../Functional Categorization/20190808_NeuronFunctionalCategorization_RCluster.xlsx", sheet = "Sheet1")
      
      
      if (length(indGS_inData)>1){
        
        ## C2 v AD1 ##
        pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD1_rev/geneSetHeatMap_C2vAD1_",indGeneSetPlot,".png",sep='')
        png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #png(pngPath, height=700, width=650, pointsize=12)
        print({
          
          blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
          rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
          blah=blah[!is.na(blah[,1]),]
          
          indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
          
          rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
          
          
          hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
          hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
          
          #RedGreenDiff=rowMeans(blah[,indRed])-rowMeans(blah[,indGreen])
          RedGreenDiff=rowMeans(blah[,indGreen])-rowMeans(blah[,indRed])
          
          
          indBlah.RG=sort(RedGreenDiff, index.return=TRUE)$ix
          
          heatmap3(blah[indBlah.RG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                   col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                   Rowv=NA, Colv=NA,  scale="row",
                   labRow = rowLabels[indBlah.RG], labCol = NA, cexRow=0.15,
                   main=rownames(geneSetEnrich)[indGeneSetPlot]) 
        })
        dev.off()
        
        ## C2 v AD2 ##
        pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD2_rev/geneSetHeatMap_C2vAD2_",indGeneSetPlot,".png",sep='')
        png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #png(pngPath, height=700, width=650, pointsize=12)
        print({
          
          blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
          rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
          blah=blah[!is.na(blah[,1]),]
          
          indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
          
          rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
          
          
          hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
          hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
          
          PurpleGreenDiff=rowMeans(blah[,indGreen])-rowMeans(blah[,indPurple])
          
          indBlah.PG=sort(PurpleGreenDiff, index.return=TRUE)$ix
          
          heatmap3(blah[indBlah.PG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                   col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                   Rowv=NA, Colv=NA,  scale="row",
                   labRow = rowLabels[indBlah.PG], labCol = NA, cexRow=0.15,
                   main=rownames(geneSetEnrich)[indGeneSetPlot]) 
        })
        dev.off()
        
        ## C2 v MCI1 ##
        pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI1_rev/geneSetHeatMap_C2vMCI1_",indGeneSetPlot,".png",sep='')
        png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #png(pngPath, height=700, width=650, pointsize=12)
        print({
          
          blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
          rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
          blah=blah[!is.na(blah[,1]),]
          
          indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
          
          rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
          
          
          hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
          hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
          
          PapayaGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indPlum])
          
          indBlah.PapG =sort(PapayaGreenDiff, index.return=TRUE)$ix
          
          heatmap3(blah[indBlah.PapG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                   col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                   Rowv=NA, Colv=NA,  scale="row",
                   labRow = rowLabels[indBlah.PapG], labCol = NA, cexRow=0.15,
                   main=rownames(geneSetEnrich)[indGeneSetPlot]) 
        })
        dev.off()
        
        ## C2 v MCI2 ##
        pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI2_rev/geneSetHeatMap_C2vMCI2_",indGeneSetPlot,".png",sep='')
        png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #png(pngPath, height=700, width=650, pointsize=12)
        print({
          
          blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
          rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
          blah=blah[!is.na(blah[,1]),]
          
          indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
          
          rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
          
          
          hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
          hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
          
          IndianRedGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indIndianred])
          
          indBlah.IRG =sort(IndianRedGreenDiff, index.return=TRUE)$ix
          
          heatmap3(blah[indBlah.IRG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                   col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                   Rowv=NA, Colv=NA,  scale="row",
                   labRow = rowLabels[indBlah.IRG], labCol = NA, cexRow=0.15,
                   main=rownames(geneSetEnrich)[indGeneSetPlot]) 
        })
        dev.off()
        
        ## C2 v C1 ##
        pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".png",sep='')
        png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #png(pngPath, height=700, width=650, pointsize=12)
        print({
          
          blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
          rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
          blah=blah[!is.na(blah[,1]),]
          
          indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
          
          rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
          
          
          hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
          hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
          
          BlueGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indBlue])
          
          indBlah.BG =sort(BlueGreenDiff, index.return=TRUE)$ix
          
          heatmap3(blah[indBlah.BG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                   col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                   Rowv=NA, Colv=NA,  scale="row",
                   labRow = rowLabels[indBlah.BG], labCol = NA, cexRow=0.15,
                   main=rownames(geneSetEnrich)[indGeneSetPlot]) 
        })
        dev.off()
        
        ## MCI2 v MCI1 ##
        pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_MCI2vMCI1_rev/geneSetHeatMap_MCI2vMCI1_",indGeneSetPlot,".png",sep='')
        png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #png(pngPath, height=700, width=650, pointsize=12)
        print({
          
          blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
          rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
          blah=blah[!is.na(blah[,1]),]
          
          indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
          
          rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
          
          
          hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
          hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
          
          IRPlumDiff= rowMeans(blah[,indIndianred])-rowMeans(blah[,indPlum])
          
          indBlah.IRP =sort(IRPlumDiff, index.return=TRUE)$ix
          
          heatmap3(blah[indBlah.IRP,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                   col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                   Rowv=NA, Colv=NA,  scale="row",
                   labRow = rowLabels[indBlah.IRP], labCol = NA, cexRow=0.15,
                   main=rownames(geneSetEnrich)[indGeneSetPlot]) 
        })
        dev.off()
        
      
        
      }
      graphics.off()
      
      RedGreenDiffOut=cbind(RedGreenDiff[indBlah.RG[length(indBlah.RG):1]],rowLabels[indBlah.RG[length(indBlah.RG):1]])
      PurpleGreenDiffOut=cbind(PurpleGreenDiff[indBlah.PG[length(indBlah.PG):1]],rowLabels[indBlah.PG[length(indBlah.PG):1]])
      PapayaGreenDiffOut=cbind(PapayaGreenDiff[indBlah.PapG[length(indBlah.PapG):1]],rowLabels[indBlah.PapG[length(indBlah.PapG):1]])
      IndianRedGreenDiffOut=cbind(IndianRedGreenDiff[indBlah.IRG[length(indBlah.IRG):1]],rowLabels[indBlah.IRG[length(indBlah.IRG):1]])
      BlueGreenDiffOut=cbind(BlueGreenDiff[indBlah.BG[length(indBlah.BG):1]],rowLabels[indBlah.BG[length(indBlah.BG):1]])
      IRPlumDiffOut=cbind(IRPlumDiff[indBlah.IRP[length(indBlah.IRP):1]],rowLabels[indBlah.IRP[length(indBlah.IRP):1]])
      
      write.xlsx(RedGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD1_rev/geneSetHeatMap_C2vAD1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
      write.xlsx(PurpleGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD2_rev/geneSetHeatMap_C2vAD2_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
      write.xlsx(PapayaGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI1_rev/geneSetHeatMap_C2vMCI1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
      write.xlsx(IndianRedGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI2_rev/geneSetHeatMap_C2vMCI2_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
      write.xlsx(BlueGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
      write.xlsx(IRPlumDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
      
    }
    
    
    #Make ES difference bar plot for neurons
    
    
    ### C2 v AD1
    geneSetDifferenceData.C2vA1=data.frame(meanBluevsOR,rownames(geneSetEnrich),pBootFDRNeu[,4])
    
    colnames(geneSetDifferenceData.C2vA1)=c("difference","geneSetName", 'pBootFDR')
    
    #Fix the spelling of "miscellaneous
    geneSetDifferenceData.C2vA1[which(geneSetDifferenceData.C2vA1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
    
   
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons_C2vAD1_rev/ES_diff_FXN_Neu_C2vAD1_",geneSetNumber,".pdf",sep='')
    #png(pngPath, width=6,height=4,units="in",res=600)
    pdf(pngPath, width=6,height=4)
    
    print({
      ggplot(geneSetDifferenceData.C2vA1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                              fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
        
        geom_bar(position=position_dodge(), stat="identity", colour="black") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
              text = element_text(size=20),
              axis.text.x = element_text(face = "plain", color = "black", 
                                         size = 18),
              axis.text.y = element_text(face = "plain", color = "black", 
                                         size = 16),
              legend.position = "none")+
        scale_fill_manual(name = "area", values=c("red","grey50"))+
        xlab("Gene Set")+
        ylab("ES Diff. (C2-A1)")+
        geom_hline(yintercept=0)+
        geom_text(aes(label=format((pBootFDRNeu[,4]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
        geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
        coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vA1)[1]+2))
      
    })
    dev.off()
    
    ### C2 v AD2
    geneSetDifferenceData.C2vA2=data.frame(meanBluevsOR4,rownames(geneSetEnrich),pBootFDRNeu[,5])
    
    colnames(geneSetDifferenceData.C2vA2)=c("difference","geneSetName", 'pBootFDR')
    
    #Fix the spelling of "miscellaneous
    geneSetDifferenceData.C2vA2[which(geneSetDifferenceData.C2vA2$geneSetName=="Miscelaneous"),2]="Miscellaneous"
    
    
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons_C2vAD2_rev/ES_diff_FXN_Neu_C2vAD2_",geneSetNumber,".pdf",sep='')
    #png(pngPath, width=6,height=4,units="in",res=600)
    pdf(pngPath, width=6,height=4)
    
    print({
      ggplot(geneSetDifferenceData.C2vA2, aes(x=reorder(geneSetName,-difference), y=difference, 
                                              fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
        
        geom_bar(position=position_dodge(), stat="identity", colour="black") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
              text = element_text(size=20),
              axis.text.x = element_text(face = "plain", color = "black", 
                                         size = 18),
              axis.text.y = element_text(face = "plain", color = "black", 
                                         size = 16),
              legend.position = "none")+
        scale_fill_manual(name = "area", values=c("red","grey50"))+
        xlab("Gene Set")+
        ylab("ES Diff. (C2-A2)")+
        geom_hline(yintercept=0)+
        geom_text(aes(label=format((pBootFDRNeu[,5]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
        geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
        coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vA2)[1]+2))
      
    })
    dev.off()
    
    ### C2 v MCI1
    geneSetDifferenceData.C2vMCI1=data.frame(meanBluevsPapaya,rownames(geneSetEnrich),pBootFDRNeu[,2])
    
    colnames(geneSetDifferenceData.C2vMCI1)=c("difference","geneSetName", 'pBootFDR')
    
    #Fix the spelling of "miscellaneous
    geneSetDifferenceData.C2vMCI1[which(geneSetDifferenceData.C2vMCI1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
    
    
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons_C2vMCI1_rev/ES_diff_FXN_Neu_C2vMCI1_",geneSetNumber,".pdf",sep='')
    #png(pngPath, width=6,height=4,units="in",res=600)
    pdf(pngPath, width=6,height=4)
    
    print({
      ggplot(geneSetDifferenceData.C2vMCI1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                              fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
        
        geom_bar(position=position_dodge(), stat="identity", colour="black") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
              text = element_text(size=20),
              axis.text.x = element_text(face = "plain", color = "black", 
                                         size = 18),
              axis.text.y = element_text(face = "plain", color = "black", 
                                         size = 16),
              legend.position = "none")+
        scale_fill_manual(name = "area", values=c("red","grey50"))+
        xlab("Gene Set")+
        ylab("ES Diff. (C2-MCI1)")+
        geom_hline(yintercept=0)+
        geom_text(aes(label=format((pBootFDRNeu[,2]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
        geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
        coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vMCI1)[1]+2))
      
    })
    dev.off()
    
    ### C2 v MCI2
    geneSetDifferenceData.C2vMCI2=data.frame(meanBluevsPeach,rownames(geneSetEnrich),pBootFDRNeu[,3])
    
    colnames(geneSetDifferenceData.C2vMCI2)=c("difference","geneSetName", 'pBootFDR')
    
    #Fix the spelling of "miscellaneous
    geneSetDifferenceData.C2vMCI2[which(geneSetDifferenceData.C2vMCI2$geneSetName=="Miscelaneous"),2]="Miscellaneous"
    
   
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons_C2vMCI2_rev/ES_diff_FXN_Neu_C2vMCI2_",geneSetNumber,".pdf",sep='')
    #png(pngPath, width=6,height=4,units="in",res=600)
    pdf(pngPath, width=6,height=4)
    
    print({
      ggplot(geneSetDifferenceData.C2vMCI2, aes(x=reorder(geneSetName,-difference), y=difference, 
                                              fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
        
        geom_bar(position=position_dodge(), stat="identity", colour="black") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
              text = element_text(size=20),
              axis.text.x = element_text(face = "plain", color = "black", 
                                         size = 18),
              axis.text.y = element_text(face = "plain", color = "black", 
                                         size = 16),
              legend.position = "none")+
        #scale_fill_manual(name = "area", values=c("red","grey50"))+
        scale_fill_manual(name = "area", values=c("grey50"))+
        xlab("Gene Set")+
        ylab("ES Diff. (C2-MCI2)")+
        geom_hline(yintercept=0)+
        geom_text(aes(label=format((pBootFDRNeu[,3]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
        geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
        coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vMCI2)[1]+2))
      
    })
    dev.off()
    
    ### C2 v C1
    geneSetDifferenceData.C2vC1=data.frame(meanBluevsBlue4,rownames(geneSetEnrich),pBootFDRNeu[,1])
    
    colnames(geneSetDifferenceData.C2vC1)=c("difference","geneSetName", 'pBootFDR')
    
    #Fix the spelling of "miscellaneous
    geneSetDifferenceData.C2vC1[which(geneSetDifferenceData.C2vC1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
    

    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons_C2vC1_rev/ES_diff_FXN_Neu_C2vC1_",geneSetNumber,".pdf",sep='')
    #png(pngPath, width=6,height=4,units="in",res=600)
    pdf(pngPath, width=6,height=4)
    
    print({
      ggplot(geneSetDifferenceData.C2vC1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                              fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
        
        geom_bar(position=position_dodge(), stat="identity", colour="black") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
              text = element_text(size=20),
              axis.text.x = element_text(face = "plain", color = "black", 
                                         size = 18),
              axis.text.y = element_text(face = "plain", color = "black", 
                                         size = 16),
              legend.position = "none")+
        scale_fill_manual(name = "area", values=c("red","grey50"))+
        xlab("Gene Set")+
        ylab("ES Diff. (C2-C1)")+
        geom_hline(yintercept=0)+
        geom_text(aes(label=format((pBootFDRNeu[,1]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
        geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
        coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]+2))
      
    })
    dev.off()
    
    ### MCI2 v MCI1
    geneSetDifferenceData.MCI2vMCI1=data.frame(meanPeachvsPapaya,rownames(geneSetEnrich),pBootFDRNeu[,6])
    
    colnames(geneSetDifferenceData.MCI2vMCI1)=c("difference","geneSetName", 'pBootFDR')
    
    #Fix the spelling of "miscellaneous
    geneSetDifferenceData.MCI2vMCI1[which(geneSetDifferenceData.MCI2vMCI1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
    
    
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Neurons_MCI2vMCI1_rev/ES_diff_FXN_Neu_MCI2vMCI1_",geneSetNumber,".pdf",sep='')
    #png(pngPath, width=6,height=4,units="in",res=600)
    pdf(pngPath, width=6,height=4)
    
    print({
      ggplot(geneSetDifferenceData.MCI2vMCI1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                              fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
        
        geom_bar(position=position_dodge(), stat="identity", colour="black") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
              text = element_text(size=20),
              axis.text.x = element_text(face = "plain", color = "black", 
                                         size = 18),
              axis.text.y = element_text(face = "plain", color = "black", 
                                         size = 16),
              legend.position = "none")+
        scale_fill_manual(name = "area", values=c("red","grey50"))+
        xlab("Gene Set")+
        ylab("ES Diff. (MCI2-MCI1)")+
        geom_hline(yintercept=0)+
        geom_text(aes(label=format((pBootFDRNeu[,6]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
        geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
        coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]+2))
      
    })
    dev.off()   
    
    
    
    
    
#____Astrocytes####
    geneSetsMat=read_excel(paste("../Functional Categorization/GeneSets_Astrocytes_FXN", geneSetNumber, ".xlsx",sep=''), sheet="sheet1")
    geneSets=as.list(geneSetsMat)
    geneSets=lapply(geneSets, function(x) x[!is.na(x)])
    
    geneSetEnrich=gsva(yDataforGSVA[clusterNumberCTDClean==4,], geneSets, mx.diff=TRUE, kcdf="Gaussian")
    
    geneSetEnrich_Astrocytes=geneSetEnrich
    
    indControl <-  sort(c(indBlue, indGreen))
    indMCI <- sort(c(indPlum, indIndianred))
    indAD <- sort(c(indRed, indPurple))
    
    #----GSVA vs Pathology----
    
    #Braak
    pathModel_pVal=matrix(NA,nrow=dim(geneSetEnrich)[1], ncol = 1)
    
      for (i in 1:dim(geneSetEnrich)[1]){
      
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,indIndianred])
      pathValues=as.numeric(pathScores[1,indIndianred])
      pathModel=lm(GSVAScores ~ pathValues)
      pathModel_pVal[i]=summary(pathModel)$coefficients[2,4]
      
      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
      
      pngPath=paste("./Braak_Astrocyte/MCI2/GSVA_vs_Braak_",i,".png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({
        
        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1) +
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "bold", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "bold", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("Braak Score")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])+
          geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])+
          annotate("text", x = -Inf, y = -Inf, label="test" , hjust = -0.5, vjust = -1)
      })
      dev.off()
      
        }
    pathModelOut=data.frame(geneSet=rownames(geneSetEnrich),pValSlope=pathModel_pVal)
    write.xlsx(pathModelOut, file="BraakvsAstGeneSets_MCI2.xlsx", sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)
    
    
    #CERAD
    pathModel_pVal=matrix(NA,nrow=dim(geneSetEnrich)[1], ncol = 1)
    
    for (i in 1:dim(geneSetEnrich)[1]){
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,indPurple])
      pathValues=as.numeric(pathScores[2,indPurple])
      pathModel=lm(GSVAScores ~ pathValues)
      pathModel_pVal[i]=summary(pathModel)$coefficients[2,4]
      
      
      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
      
      pngPath=paste("./CERAD_Astrocyte/AD2/GSVA_vs_CERAD_",i,".png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({
        
        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1) +
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "plain", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("CERAD")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])+
          geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
      })
      dev.off()
      
    }
    pathModelOut=data.frame(geneSet=rownames(geneSetEnrich),pValSlope=pathModel_pVal)
    write.xlsx(pathModelOut, file="CERADvsAstGeneSets_AD2.xlsx", sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)
    
    
    #APOE
    for (i in 1:dim(geneSetEnrich)[1]){
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,])
      pathValues=as.numeric(pathScores[3,])
      
      GSVAScores=GSVAScores[-which(is.na(pathValues))]
      pathValues=pathValues[-which(is.na(pathValues))]
      
      pathModel=lm(GSVAScores ~ pathValues)
      
      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
      
      pngPath=paste("./FXN1_norm_Astrocytes/GSVA_vs_APOE_",i,".png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({
        
        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1) +
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "plain", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("APOE")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])

      })
      dev.off()
      
    }
    
    #Sex
    for (i in 1:dim(geneSetEnrich)[1]){
      #Format data for regression
      GSVAScores=as.numeric(geneSetEnrich[i,])
      pathValues=as.numeric(pathScores[4,])
      pathModel=lm(GSVAScores ~ pathValues)
      
      dataPath=data.frame(pathValues,GSVAScores)
      statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
      
      pngPath=paste("./FXN1_norm_Astrocytes/GSVA_vs_Sex_",i,".png",sep='')
      png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
      print({
        
        e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
        e+ geom_violin()+
          geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                        width=0.4, colour="black", alpha=1, size=1) +
          geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                        size=0.5,col="red", width = .5)+
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black",
                                           size = 16),
                axis.text.y = element_text(face = "plain", color = "black",
                                           size = 16),
                legend.position = "none", plot.title = element_text(hjust = 0.5))+
          xlab("Sex")+
          ylab("Enrichment Score (a.u.)")+
          ggtitle(row.names(geneSetEnrich)[i])+
          scale_x_discrete(labels=c("0" = "Female", "1" = "Male"))+
          geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
      })
      dev.off()
      
    }
    
    #----Astrocyte clustering analysis of GSVA----
    
    geneSetEnrichSort=geneSetEnrich_Astrocytes[,indSort] #sort samples by phenotype
    geneSetEnrichSortZ=t(apply(geneSetEnrichSort,1,scale)); #z-score the data
    
    hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes/GSVA_cluster.png",sep='')
    png(pngPath, width=8,height=4,units="in",res=600, pointsize = 16)
    
    heatmap3(geneSetEnrichSortZ, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
             col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
             Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
             labRow = colnames(geneSetsMat), labCol = NA, cexRow=1.2) 
    title("GSVA for astrocyte gene sets", adj = 0.5, line = 2)
    
    dev.off()
    
    
    corGSVA=cor(t(geneSetEnrichSortZ))
    hCor=hclust(as.dist((1-corGSVA)/2))
    
    pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes/GSVA_corr.png",sep='')
    png(pngPath, width=6,height=8,units="in",res=600, pointsize = 14)
    #par(mar=c(5,6,4,4)+.1)
    
    heatmap3(corGSVA,
             col=barColorsCor, breaks=breakBarColorsCor,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
             Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",
             labRow = colnames(geneSetsMat), labCol = colnames(geneSetsMat),margins=c(12,13)) 
    title("GSVA for astrocyte gene sets", adj = 0.5, line = 1)
    dev.off()
    
    

    

    
    indclusterNumberCTDClean=which(clusterNumberCTDClean==4)
    
    blah=yDataforGSVA[indclusterNumberCTDClean,]
    
    #Permutation p-value
    if (compPermute==1){
      R=1000
      meanOut6=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut5=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut4=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut3=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut2=matrix(NA,nrow=length(geneSets), ncol = R)
      meanOut1=matrix(NA,nrow=length(geneSets), ncol = R)
      
      for (i in 1:R)
      {
        blah2=blah
        rownames(blah2)=rownames(blah)[sample(1:dim(blah)[1])]
        gsvaBoot=gsva(blah2,geneSets, mx.diff=TRUE, kcdf="Gaussian")
        meanOut6[,i]=rowMeans(gsvaBoot[,indIndianred])-rowMeans(gsvaBoot[,indPlum]) # MCI2 - MCI1
        meanOut5[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPurple]) # AD2 - C2
        meanOut4[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indRed]) # AD1 - C2 
        meanOut3[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indIndianred]) # MCI2 - C2
        meanOut2[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPlum]) #MCI1 - C2
        meanOut1[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indBlue]) # C1 - C2
        
        
        print(i)
      }
      meanTrue6=rowMeans(geneSetEnrich[,indIndianred])-rowMeans(geneSetEnrich[,indPlum])
      meanTrue5=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPurple])
      meanTrue4=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indRed])
      meanTrue3=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indIndianred])
      meanTrue2=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPlum])
      meanTrue1=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indBlue])
      
      
      pBoot6=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot5=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot4=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot3=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot2=matrix(NA,nrow=length(geneSets), ncol = 1)
      pBoot1=matrix(NA,nrow=length(geneSets), ncol = 1)
      
      for (j in 1:length(geneSets))
        
      {
        pBoot6[j] = mean(abs(meanOut6[j,]) > abs(meanTrue6[j]))
        pBoot5[j] = mean(abs(meanOut5[j,]) > abs(meanTrue5[j]))
        pBoot4[j] = mean(abs(meanOut4[j,]) > abs(meanTrue4[j]))
        pBoot3[j] = mean(abs(meanOut3[j,]) > abs(meanTrue3[j]))
        pBoot2[j] = mean(abs(meanOut2[j,]) > abs(meanTrue2[j]))
        pBoot1[j] = mean(abs(meanOut1[j,]) > abs(meanTrue1[j]))
        
      }
      pBootFDR6=p.adjust(pBoot6,method = "fdr",)
      pBootFDR5=p.adjust(pBoot5,method = "fdr",)
      pBootFDR4=p.adjust(pBoot4,method = "fdr",)
      pBootFDR3=p.adjust(pBoot3,method = "fdr",)
      pBootFDR2=p.adjust(pBoot2,method = "fdr",)
      pBootFDR1=p.adjust(pBoot1,method = "fdr",)
      
      
      pBoot=data.frame(pBoot1, pBoot2, pBoot3, pBoot4, pBoot5, pBoot6)
      pBootFDR=data.frame(pBootFDR1, pBootFDR2, pBootFDR3, pBootFDR4, pBootFDR5, pBootFDR6)
      
      rownames(pBoot)=rownames(geneSetEnrich)
      rownames(pBootFDR)=rownames(geneSetEnrich)
      #End permutation
      
      pBootAst=pBoot
      pBootFDRAst=pBootFDR
      
    }
    
    if(geneSetNumber=="1"){
      save_pBootAst1=data.frame(pBootAst,pBootFDRAst)
    }else{
      save_pBootAst2=data.frame(pBootAst,pBootFDRAst)
    }
    
    
    
    if(compPermute=="0" & geneSetNumber=="2"){
    
      pBootAst=save_pBootAst2
      }else if(compPermute=="0" & geneSetNumber=="1"){
        pBootAst=save_pBootAst1
      
    }
    
    
    indKeepGeneSets=which(lengths(geneSets, use.names = FALSE)>9)
    
    geneSetEnrichSort=geneSetEnrichSort[indKeepGeneSets,]
    geneSetEnrich=geneSetEnrich[indKeepGeneSets,]
    pBootAst=pBootAst[indKeepGeneSets,]
    pBootFDRAst = pBootFDRAst[indKeepGeneSets,]
    
    
    #plotting astrocyte genes
    heatmap3(yDataforGSVA[clusterNumberCTDClean==4,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
             col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
             Rowv=NA, Colv=NA,  scale="row",
             labRow = NA, labCol = NA, main = "Astrocyte Genes") 
    
    #Compute p-values and FDR q-values for each gene set.
    # pValsTukey=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
    # pANOVA=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
    numGenes=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
    
    for (i in 1:dim(geneSetEnrich)[1])
    { #p=t.test(geneSetEnrichSort[i,allColorGroups=="Blue"], geneSetEnrichSort[i,allColorGroups=="Purple"])
      #pVals[i]=p$p.value
      
      indGetGS=which(match(names(geneSets),row.names(geneSetEnrich)[i])==1)
      numGenes[i]=length(unlist(geneSets[indGetGS]))
      
    #   #ANOVA
    #   dataANOVA=data.frame(geneSetEnrichSort[i,], allColorGroups)
    #   colnames(dataANOVA)=c("GSVA_ES","color")
    #   AOVObject=aov(GSVA_ES ~ color, data = dataANOVA)
    #   resultANOVA=summary(AOVObject)
    #   pANOVA[i]=resultANOVA[[1]][["Pr(>F)"]][1]
    #   
    #   ANOVATukey=TukeyHSD(AOVObject,conf.level = 0.95)
    #   pValsTukey[i]=ANOVATukey$color[which(rownames(ANOVATukey$color)=="skyblue4-orangered4"),4] #Compare blue vs purple
     }
    
    # FDRANOVA=p.adjust(pANOVA,method = "fdr")
    # qFDRTukey=p.adjust(pValsTukey,method = "fdr")
    geneSets_Stats=data.frame(rownames(geneSetEnrichSort),pBootAst,pBootFDRAst, numGenes);
    write.xlsx(geneSets_Stats, file=paste("geneSetStats_Astrocytes_FXN_reduceGeneSets_C2vALL_",geneSetNumber,"_norm.xlsx",sep=''), sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)
    
    
      
      #Plot bar graphs and gene set heatmaps
    
      meanBluevsOR=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) #C2 v AD1
      meanBluevsOR4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v AD2
      meanBluevsPapaya=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v MCI1
      meanBluevsPeach=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v MCI2
      meanBluevsBlue4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v C1
      meanPeachvsPapaya =matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) #MCI2 v MCI1
      
      for (indGeneSetPlot in 1:dim(geneSetEnrich)[1])
      {
        dataBar=data.frame(as.data.frame(geneSetEnrich[indGeneSetPlot,indSort]), as.data.frame(allColorGroups))
        colnames(dataBar)=c("ES","colors")
        
        #Summarize stats for each group color
        statsOut = summarySE(dataBar, measurevar="ES", groupvars="colors")
        
        meanBluevsOR[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered"] # AD1 - C2
        meanBluevsOR4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered4"] #AD2 - C2
        meanBluevsPapaya[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="papayawhip"] #C2 v MC1
        meanBluevsPeach[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="peachpuff3"] #C2 vs MC2
        meanBluevsBlue4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="skyblue4"] # C2 vs C1
        meanPeachvsPapaya[indGeneSetPlot]=statsOut$ES[statsOut$color=="peachpuff3"]-statsOut$ES[statsOut$color=="papayawhip"] # MCI2 v MCI1
        
        #Re-order the colors to match the heatmaps
        statsOut$colors = factor(statsOut$colors, levels = c("skyblue4","skyblue","papayawhip","peachpuff3","orangered","orangered4"))
        
        group.colors <- c(skyblue4 = "skyblue4", skyblue = "skyblue", papayawhip="papayawhip",peachpuff3="peachpuff3", orangered ="orangered", orangered4 = "orangered4")
        
        pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes/geneSets",indGeneSetPlot,".png",sep='')
        png(pngPath, width=5,height=4,units="in",res=600)
        print({
          ggplot(statsOut, aes(x=colors, y=ES, fill=colors)) + 
            geom_bar(position=position_dodge(), stat="identity", colour="black") +
            geom_errorbar(aes(ymin=ES-se, ymax=ES+se),
                          width=.2, size=1.5,                    # Width of the error bars
                          position=position_dodge(.9))+
            theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                  text = element_text(size=28),
                  axis.text.x = element_text(face = "plain", color = "black", 
                                             size = 24),
                  axis.text.y = element_text(face = "plain", color = "black", 
                                             size = 24),
                  legend.position = "none")+
            xlab("")+
            ylab("Enrichment Score")+
            ggtitle(rownames(geneSetEnrich)[indGeneSetPlot])+
            scale_fill_manual(values=group.colors)+
            ylim(-0.5,0.5)+
            scale_x_discrete(labels=c("C1", "C2", "M1", "M2", "A1", "A2"))+
            geom_hline(yintercept=0)
          
        })
        dev.off()
        
        #plot heatmaps for each gene set
        
        #Get index of gene set associated with enrichment
        indGetGS=which(match(names(geneSets),row.names(geneSetEnrich)[indGeneSetPlot])==1)
        
        
        indGS_inData=match(as.vector(unlist(geneSets[indGetGS])),geneSortClean)
        indGS_inData=indGS_inData[!is.na(indGS_inData)]
        
        FXN_Annotation_Astrocytes=read_excel("../Functional Categorization/20190804_AstrocyteFunctionalCategorization_RCluster.xlsx", sheet = "Sheet1")
        
        if (length(indGS_inData)>1){
          
          ## C2 v AD1 ##
          pngPath=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vAD1_rev/geneSetHeatMap_C2vAD1_",indGeneSetPlot,".png",sep='')
          png(pngPath,width=3.25,height=3.25,units="in",res=600)
          #png(pngPath, height=700, width=650, pointsize=12)
          print({
            
            blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
            rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
            blah=blah[!is.na(blah[,1]),]
            
            indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Astrocytes$GeneSymbol)
            
            rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Astrocytes$FXN2[indFXN2_in_GS], sep = "")
            
            
            hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
            hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
            
            #RedGreenDiff=rowMeans(blah[,indRed])-rowMeans(blah[,indGreen])
            RedGreenDiff=rowMeans(blah[,indGreen])-rowMeans(blah[,indRed])
            
            
            indBlah.RG=sort(RedGreenDiff, index.return=TRUE)$ix
            
            heatmap3(blah[indBlah.RG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                     col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                     Rowv=NA, Colv=NA,  scale="row",
                     labRow = rowLabels[indBlah.RG], labCol = NA, cexRow=0.15,
                     main=rownames(geneSetEnrich)[indGeneSetPlot]) 
          })
          dev.off()
          
          ## C2 v AD2 ##
          pngPath=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vAD2_rev/geneSetHeatMap_C2vAD2_",indGeneSetPlot,".png",sep='')
          png(pngPath,width=3.25,height=3.25,units="in",res=600)
          #png(pngPath, height=700, width=650, pointsize=12)
          print({
            
            blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
            rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
            blah=blah[!is.na(blah[,1]),]
            
            indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Astrocytes$GeneSymbol)
            
            rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Astrocytes$FXN2[indFXN2_in_GS], sep = "")
            
            
            hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
            hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
            
            PurpleGreenDiff=rowMeans(blah[,indGreen])-rowMeans(blah[,indPurple])
            
            indBlah.PG=sort(PurpleGreenDiff, index.return=TRUE)$ix
            
            heatmap3(blah[indBlah.PG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                     col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                     Rowv=NA, Colv=NA,  scale="row",
                     labRow = rowLabels[indBlah.PG], labCol = NA, cexRow=0.15,
                     main=rownames(geneSetEnrich)[indGeneSetPlot]) 
          })
          dev.off()
          
          ## C2 v MCI1 ##
          pngPath=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vMCI1_rev/geneSetHeatMap_C2vMCI1_",indGeneSetPlot,".png",sep='')
          png(pngPath,width=3.25,height=3.25,units="in",res=600)
          #png(pngPath, height=700, width=650, pointsize=12)
          print({
            
            blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
            rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
            blah=blah[!is.na(blah[,1]),]
            
            indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Astrocytes$GeneSymbol)
            
            rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Astrocytes$FXN2[indFXN2_in_GS], sep = "")
            
            
            hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
            hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
            
            PapayaGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indPlum])
            
            indBlah.PapG =sort(PapayaGreenDiff, index.return=TRUE)$ix
            
            heatmap3(blah[indBlah.PapG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                     col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                     Rowv=NA, Colv=NA,  scale="row",
                     labRow = rowLabels[indBlah.PapG], labCol = NA, cexRow=0.15,
                     main=rownames(geneSetEnrich)[indGeneSetPlot]) 
          })
          dev.off()
          
          ## C2 v MCI2 ##
          pngPath=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vMCI2_rev/geneSetHeatMap_C2vMCI2_",indGeneSetPlot,".png",sep='')
          png(pngPath,width=3.25,height=3.25,units="in",res=600)
          #png(pngPath, height=700, width=650, pointsize=12)
          print({
            
            blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
            rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
            blah=blah[!is.na(blah[,1]),]
            
            indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Astrocytes$GeneSymbol)
            
            rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Astrocytes$FXN2[indFXN2_in_GS], sep = "")
            
            
            hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
            hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
            
            IndianRedGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indIndianred])
            
            indBlah.IRG =sort(IndianRedGreenDiff, index.return=TRUE)$ix
            
            heatmap3(blah[indBlah.IRG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                     col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                     Rowv=NA, Colv=NA,  scale="row",
                     labRow = rowLabels[indBlah.IRG], labCol = NA, cexRow=0.15,
                     main=rownames(geneSetEnrich)[indGeneSetPlot]) 
          })
          dev.off()
          
          ## C2 v C1 ##
          pngPath=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".png",sep='')
          png(pngPath,width=3.25,height=3.25,units="in",res=600)
          #png(pngPath, height=700, width=650, pointsize=12)
          print({
            
            blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
            rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
            blah=blah[!is.na(blah[,1]),]
            
            indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Astrocytes$GeneSymbol)
            
            rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Astrocytes$FXN2[indFXN2_in_GS], sep = "")
            
            
            hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
            hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
            
            BlueGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indBlue])
            
            indBlah.BG =sort(BlueGreenDiff, index.return=TRUE)$ix
            
            heatmap3(blah[indBlah.BG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                     col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                     Rowv=NA, Colv=NA,  scale="row",
                     labRow = rowLabels[indBlah.BG], labCol = NA, cexRow=0.15,
                     main=rownames(geneSetEnrich)[indGeneSetPlot]) 
          })
          dev.off()
          
          ## MCI2 v MCI1 ##
          pngPath=paste("./FXN", geneSetNumber,"_norm_Astrocytes_MCI2vMCI1_rev/geneSetHeatMap_MCI2vMCI1_",indGeneSetPlot,".png",sep='')
          png(pngPath,width=3.25,height=3.25,units="in",res=600)
          #png(pngPath, height=700, width=650, pointsize=12)
          print({
            
            blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale)) 
            rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
            blah=blah[!is.na(blah[,1]),]
            
            indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Astrocytes$GeneSymbol)
            
            rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Astrocytes$FXN2[indFXN2_in_GS], sep = "")
            
            
            hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
            hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
            
            IRPlumDiff= rowMeans(blah[,indIndianred])-rowMeans(blah[,indPlum])
            
            indBlah.IRP =sort(IRPlumDiff, index.return=TRUE)$ix
            
            heatmap3(blah[indBlah.IRP,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
                     col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                     Rowv=NA, Colv=NA,  scale="row",
                     labRow = rowLabels[indBlah.IRP], labCol = NA, cexRow=0.15,
                     main=rownames(geneSetEnrich)[indGeneSetPlot]) 
          })
          dev.off()
          
          
          
        }
        graphics.off()
        
        RedGreenDiffOut=cbind(RedGreenDiff[indBlah.RG[length(indBlah.RG):1]],rowLabels[indBlah.RG[length(indBlah.RG):1]])
        PurpleGreenDiffOut=cbind(PurpleGreenDiff[indBlah.PG[length(indBlah.PG):1]],rowLabels[indBlah.PG[length(indBlah.PG):1]])
        PapayaGreenDiffOut=cbind(PapayaGreenDiff[indBlah.PapG[length(indBlah.PapG):1]],rowLabels[indBlah.PapG[length(indBlah.PapG):1]])
        IndianRedGreenDiffOut=cbind(IndianRedGreenDiff[indBlah.IRG[length(indBlah.IRG):1]],rowLabels[indBlah.IRG[length(indBlah.IRG):1]])
        BlueGreenDiffOut=cbind(BlueGreenDiff[indBlah.BG[length(indBlah.BG):1]],rowLabels[indBlah.BG[length(indBlah.BG):1]])
        IRPlumDiffOut=cbind(IRPlumDiff[indBlah.IRP[length(indBlah.IRP):1]],rowLabels[indBlah.IRP[length(indBlah.IRP):1]])
        
        write.xlsx(RedGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vAD1_rev/geneSetHeatMap_C2vAD1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        write.xlsx(PurpleGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vAD2_rev/geneSetHeatMap_C2vAD2_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        write.xlsx(PapayaGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vMCI1_rev/geneSetHeatMap_C2vMCI1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        write.xlsx(IndianRedGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vMCI2_rev/geneSetHeatMap_C2vMCI2_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        write.xlsx(BlueGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        write.xlsx(IRPlumDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Astrocytes_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        
      }
    
      #Make ES difference bar plot for neurons
      
      
      ### C2 v AD1
      geneSetDifferenceData.C2vA1=data.frame(meanBluevsOR,rownames(geneSetEnrich),pBootFDRAst[,4])
      
      colnames(geneSetDifferenceData.C2vA1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      geneSetDifferenceData.C2vA1[which(geneSetDifferenceData.C2vA1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes_C2vAD1_rev/ES_diff_FXN_Ast_C2vAD1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vA1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-A1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRAst[,4]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vA1)[1]+2))
        
      })
      dev.off()
      
      ### C2 v AD2
      geneSetDifferenceData.C2vA2=data.frame(meanBluevsOR4,rownames(geneSetEnrich),pBootFDRAst[,5])
      
      colnames(geneSetDifferenceData.C2vA2)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      geneSetDifferenceData.C2vA2[which(geneSetDifferenceData.C2vA2$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes_C2vAD2_rev/ES_diff_FXN_Ast_C2vAD2_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vA2, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-A2)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRAst[,5]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vA2)[1]+2))
        
      })
       dev.off()
      
      ### C2 v MCI1
      geneSetDifferenceData.C2vMCI1=data.frame(meanBluevsPapaya,rownames(geneSetEnrich),pBootFDRAst[,2])
      
      colnames(geneSetDifferenceData.C2vMCI1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      geneSetDifferenceData.C2vMCI1[which(geneSetDifferenceData.C2vMCI1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes_C2vMCI1_rev/ES_diff_FXN_Neu_C2vMCI1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vMCI1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                  fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-MCI1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRAst[,2]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vMCI1)[1]+2))
        
      })
      dev.off()
      
      ### C2 v MCI2
      geneSetDifferenceData.C2vMCI2=data.frame(meanBluevsPeach,rownames(geneSetEnrich),pBootFDRAst[,3])
      
      colnames(geneSetDifferenceData.C2vMCI2)=c("difference","geneSetName", 'pBootFDR')
  
      #Fix the spelling of "miscellaneous
      geneSetDifferenceData.C2vMCI2[which(geneSetDifferenceData.C2vMCI2$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes_C2vMCI2_rev/ES_diff_FXN_Ast_C2vMCI2_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vMCI2, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                  fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          #scale_fill_manual(name = "area", values=c("red","grey50"))+
          scale_fill_manual(name = "area", values=c("grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-MCI2)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRAst[,3]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vMCI2)[1]+2))
        
      })
      dev.off()
      
      ### C2 v C1
      geneSetDifferenceData.C2vC1=data.frame(meanBluevsBlue4,rownames(geneSetEnrich),pBootFDRAst[,1])
      
      colnames(geneSetDifferenceData.C2vC1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      geneSetDifferenceData.C2vC1[which(geneSetDifferenceData.C2vC1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes_C2vC1_rev/ES_diff_FXN_Ast_C2vC1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vC1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-C1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRAst[,1]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]+2))
        
      })
      dev.off()
      
      ### MCI2 v MCI1
      geneSetDifferenceData.MCI2vMCI1=data.frame(meanPeachvsPapaya,rownames(geneSetEnrich),pBootFDRAst[,6])
      
      colnames(geneSetDifferenceData.MCI2vMCI1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      geneSetDifferenceData.MCI2vMCI1[which(geneSetDifferenceData.MCI2vMCI1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Astrocytes_MCI2vMCI1_rev/ES_diff_FXN_Ast_MCI2vMCI1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.MCI2vMCI1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                    fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (MCI2-MCI1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRAst[,6]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]+2))
        
      })
      dev.off()   
    
    
    
    
    
    
    
      ####____MICROGLIA_____####
      
      # Extract the gene names of microglia dataset for each identified subgroup 

      
      microDataClust.1 <- ySortZMicroglia_ann[clust.microglia.cut == 1,]
      microDataClust.2 <- ySortZMicroglia_ann[clust.microglia.cut == 2,]
      microDataClust.3 <- ySortZMicroglia_ann[clust.microglia.cut == 3,]
      microDataClust.4 <- ySortZMicroglia_ann[clust.microglia.cut == 4,]
      microDataClust.5 <- ySortZMicroglia_ann[clust.microglia.cut == 5,]
      
      writeMicroData <- as.data.frame(ySortZMicroglia_ann)
      writeMicroData$symbol <- rownames(writeMicroData)
      writeMicroData$cluster <- "na"
      
      writeMicroData[clust.microglia.cut == 1,]$cluster <- "cluster1"
      writeMicroData[clust.microglia.cut == 2,]$cluster <- "cluster2"
      writeMicroData[clust.microglia.cut == 3,]$cluster <- "cluster3"
      writeMicroData[clust.microglia.cut == 4,]$cluster <- "cluster4"
      writeMicroData[clust.microglia.cut == 5,]$cluster <- "cluster5"
      
      writeMicroData <- (writeMicroData[,c("symbol","cluster")])
      
      write.xlsx(writeMicroData, file=paste("./MicrogliaClusterGeneSet.xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)
      
      
      microGenesClust.1 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 1,])
      microGenesClust.2 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 2,])
      microGenesClust.3 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 3,])
      microGenesClust.4 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 4,])
      microGenesClust.5 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 5,])
      
      
      
      ### GO analysis to identify Biological Functions of the clusters ###
      
      #convert gene symbols to EntrezIDs 
      entrezClust1 <- mapIds(org.Hs.eg.db, microGenesClust.1, "ENTREZID", "SYMBOL")
      entrezClust2 <- mapIds(org.Hs.eg.db, microGenesClust.2, "ENTREZID", "SYMBOL")
      entrezClust3 <- mapIds(org.Hs.eg.db, microGenesClust.3, "ENTREZID", "SYMBOL")
      entrezClust4 <- mapIds(org.Hs.eg.db, microGenesClust.4, "ENTREZID", "SYMBOL")
      entrezClust5 <- mapIds(org.Hs.eg.db, microGenesClust.5, "ENTREZID", "SYMBOL")
      
      #perfrom GO 
      
      GOclust1 <- enrichGO(entrezClust1, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE)
      GOclust2 <- enrichGO(entrezClust2, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE)
      GOclust3 <- enrichGO(entrezClust3, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE)
      GOclust4 <- enrichGO(entreGOzClust4, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE)
      GOclust5 <- enrichGO(entrezClust5, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE)
      
      # GOclust1 <- enrichPathway(entrezClust1, pvalueCutoff = 0.05, pAdjustMethod = "fdr", readable= TRUE)
      # GOclust2 <- enrichGO(entrezClust2, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE)
      # GOclust3 <- enrichGO(entrezClust3, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE)
      # GOclust4 <- enrichGO(entreGOzClust4, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE)
      # GOclust5 <- enrichGO(entrezClust5, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE)
      # 
      #generate geneSet for GSVA
      microGenes.total <- list(microGenesClust.1,microGenesClust.2,microGenesClust.3,microGenesClust.4,microGenesClust.5)
      list.names <- c("M1","M2","M3","M4","M5")
      names(microGenes.total) <- list.names
      microGeneSet =lapply(microGenes.total, function(x) x[!is.na(x)])
      
      # DATA for GSVA in Microglia 
      DataGSVA.microglia <- ySortNormMicrogliaGenes #are the normalized data (after vsn) but only of microglia (clust = 17) // NOT z-scored
      rownames(DataGSVA.microglia)  <- microGeneNames
      
      
      geneSetEnrich=gsva(DataGSVA.microglia, microGeneSet, mx.diff=TRUE, kcdf="Gaussian")
      
      geneSetEnrich_Microglia=geneSetEnrich
      
      #----GSVA vs Pathology----
      
      #Braak
      
      indControl <-  sort(c(indBlue, indGreen))
      indMCI <- sort(c(indPlum, indIndianred))
      indAD <- sort(c(indRed, indPurple))
      
      for (i in 1:dim(geneSetEnrich)[1]){
        
        #Format data for regression
        GSVAScores=as.numeric(geneSetEnrich[i,indGreen])
        pathValues=as.numeric(pathScores[1,indGren])
        pathModel=lm(GSVAScores ~ pathValues)
        
        dataPath=data.frame(pathValues,GSVAScores)
        statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
        
        pngPath=paste("./FXN1_norm_Neurons_GSVAvBraak_AD/GSVA_vs_Braak_AD_",i,".png",sep='')
        png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
        print({
          
          e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
          e+ geom_violin()+
            geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                          width=0.4, colour="black", alpha=1, size=1) +
            geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                          size=0.5,col="red", width = .5)+
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  text = element_text(size=20),
                  axis.text.x = element_text(face = "bold", color = "black",
                                             size = 16),
                  axis.text.y = element_text(face = "bold", color = "black",
                                             size = 16),
                  legend.position = "none", plot.title = element_text(hjust = 0.5))+
            xlab("Braak Score")+
            ylab("Enrichment Score (a.u.)")+
            ggtitle(row.names(geneSetEnrich)[i])+
            geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
        })
        dev.off()
        
      }
      
      #CERAD
      for (i in 1:dim(geneSetEnrich)[1]){
        #Format data for regression
        GSVAScores=as.numeric(geneSetEnrich[i,indAD])
        pathValues=as.numeric(pathScores[2,])
        pathModel=lm(GSVAScores ~ pathValues)
        
        dataPath=data.frame(pathValues,GSVAScores)
        statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
        
        pngPath=paste("./FXN1_norm_Neurons_GSVAvCERAD_AD/GSVA_vs_CERAD_AD_",i,".png",sep='')
        png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
        print({
          
          e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
          e+ geom_violin()+
            geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                          width=0.4, colour="black", alpha=1, size=1) +
            geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                          size=0.5,col="red", width = .5)+
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  text = element_text(size=20),
                  axis.text.x = element_text(face = "plain", color = "black",
                                             size = 16),
                  axis.text.y = element_text(face = "plain", color = "black",
                                             size = 16),
                  legend.position = "none", plot.title = element_text(hjust = 0.5))+
            xlab("CERAD")+
            ylab("Enrichment Score (a.u.)")+
            ggtitle(row.names(geneSetEnrich)[i])+
            geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
        })
        dev.off()
        
      }
      
      #Plot CERAD vs Braak
      x=as.matrix(pathScores[1,])
      y=as.matrix(pathScores[2,])
      
      
      plot(x,y, xlab="Braak", 
           ylab="CERAD")
      
      
      #APOE
      for (i in 1:dim(geneSetEnrich)[1]){
        #Format data for regression
        GSVAScores=as.numeric(geneSetEnrich[i,])
        pathValues=as.numeric(pathScores[3,])
        
        GSVAScores=GSVAScores[-which(is.na(pathValues))]
        pathValues=pathValues[-which(is.na(pathValues))]
        
        pathModel=lm(GSVAScores ~ pathValues)
        
        dataPath=data.frame(pathValues,GSVAScores)
        statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
        
        pngPath=paste("./FXN1_norm_Neurons/GSVA_vs_APOE_",i,".png",sep='')
        png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
        print({
          
          e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
          e+ geom_violin()+
            geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                          width=0.4, colour="black", alpha=1, size=1) +
            geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                          size=0.5,col="red", width = .5)+
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  text = element_text(size=20),
                  axis.text.x = element_text(face = "plain", color = "black",
                                             size = 16),
                  axis.text.y = element_text(face = "plain", color = "black",
                                             size = 16),
                  legend.position = "none", plot.title = element_text(hjust = 0.5))+
            xlab("APOE")+
            ylab("Enrichment Score (a.u.)")+
            ggtitle(row.names(geneSetEnrich)[i])
        })
        dev.off()
        
      }
      
      #Sex
      for (i in 1:dim(geneSetEnrich)[1]){
        #Format data for regression
        GSVAScores=as.numeric(geneSetEnrich[i,])
        pathValues=as.numeric(pathScores[4,])
        pathModel=lm(GSVAScores ~ pathValues)
        
        dataPath=data.frame(pathValues,GSVAScores)
        statsOut = summarySE(dataPath, measurevar="GSVAScores", groupvars="pathValues")
        
        pngPath=paste("./FXN1_norm_Neurons/GSVA_vs_Sex_",i,".png",sep='')
        png(pngPath, width=5,height=4,units="in",res=600, pointsize = 12)
        print({
          
          e=ggplot(data=dataPath, aes(x = as.factor(pathValues), y = GSVAScores))
          e+ geom_violin()+
            geom_errorbar(data=statsOut, aes(x=as.factor(pathValues), ymin=GSVAScores-se, ymax=GSVAScores+se),
                          width=0.4, colour="black", alpha=1, size=1) +
            geom_crossbar(data=statsOut, aes(x=as.factor(pathValues),ymin = GSVAScores, ymax = GSVAScores),
                          size=0.5,col="red", width = .5)+
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  text = element_text(size=20),
                  axis.text.x = element_text(face = "plain", color = "black",
                                             size = 16),
                  axis.text.y = element_text(face = "plain", color = "black",
                                             size = 16),
                  legend.position = "none", plot.title = element_text(hjust = 0.5))+
            xlab("Sex")+
            ylab("Enrichment Score (a.u.)")+
            ggtitle(row.names(geneSetEnrich)[i])+
            scale_x_discrete(labels=c("0" = "Female", "1" = "Male"))+
            geom_abline(intercept=pathModel$coefficients[1], slope=pathModel$coefficients[2])
        })
        dev.off()
        
      }
      #----GSVA heatmap----
      
      
      geneSetEnrichSort=geneSetEnrich_Microglia[,indSort] #sort samples by phenotype
      geneSetEnrichSortZ=t(apply(geneSetEnrich_Microglia,1,scale)); #z-score the data
      
      hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia/GSVA_cluster.png",sep='')
      png(pngPath, width=8,height=4,units="in",res=600, pointsize = 16)
      
      heatmap3(geneSetEnrichSortZ, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
               col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
               Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
               labRow = colnames(microGeneSet), labCol = NA, cexRow=1.2) 
      
      title("GSVA for microglia gene sets", adj = 0.5, line = 2)
      
      dev.off()
      
      
      corGSVA=cor(t(geneSetEnrichSortZ))
      hCor=hclust(as.dist((1-corGSVA)/2))
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia/GSVA_corr.png",sep='')
      png(pngPath, width=6,height=8,units="in",res=600, pointsize = 14)
      #par(mar=c(5,6,4,4)+.1)
      
      heatmap3(corGSVA,
               col=barColorsCor, breaks=breakBarColorsCor,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
               Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",
               labRow = colnames(microGeneSet), labCol = colnames(microGeneSet),margins=c(12,13)) 
      title("GSVA for microglia gene sets", adj = 0.5, line = 1)
      dev.off()
      
      blah=DataGSVA.microglia
       
        
        #Permutation p-value
        if (compPermute==1){
          R=1000
          meanOut6=matrix(NA,nrow=length(microGeneSet), ncol = R)
          meanOut5=matrix(NA,nrow=length(microGeneSet), ncol = R)
          meanOut4=matrix(NA,nrow=length(microGeneSet), ncol = R)
          meanOut3=matrix(NA,nrow=length(microGeneSet), ncol = R)
          meanOut2=matrix(NA,nrow=length(microGeneSet), ncol = R)
          meanOut1=matrix(NA,nrow=length(microGeneSet), ncol = R)
          
          for (i in 1:R)
          {
            blah2=blah
            rownames(blah2)=rownames(blah)[sample(1:dim(blah)[1])]
            gsvaBoot=gsva(blah2,microGeneSet, mx.diff=TRUE, kcdf="Gaussian")
            meanOut6[,i]=rowMeans(gsvaBoot[,indIndianred])-rowMeans(gsvaBoot[,indPlum]) # MCI2 - MCI1
            meanOut5[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPurple]) # AD2 - C2
            meanOut4[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indRed]) # AD1 - C2 
            meanOut3[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indIndianred]) # MCI2 - C2
            meanOut2[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPlum]) #MCI1 - C2
            meanOut1[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indBlue]) # C1 - C2
            
            
            print(i)
          }
          meanTrue6=rowMeans(geneSetEnrich[,indIndianred])-rowMeans(geneSetEnrich[,indPlum])
          meanTrue5=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPurple])
          meanTrue4=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indRed])
          meanTrue3=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indIndianred])
          meanTrue2=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPlum])
          meanTrue1=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indBlue])
          
          
          pBoot6=matrix(NA,nrow=length(microGeneSet), ncol = 1)
          pBoot5=matrix(NA,nrow=length(microGeneSet), ncol = 1)
          pBoot4=matrix(NA,nrow=length(microGeneSet), ncol = 1)
          pBoot3=matrix(NA,nrow=length(microGeneSet), ncol = 1)
          pBoot2=matrix(NA,nrow=length(microGeneSet), ncol = 1)
          pBoot1=matrix(NA,nrow=length(microGeneSet), ncol = 1)
          
          for (j in 1:length(microGeneSet))
            
          {
            pBoot6[j] = mean(abs(meanOut6[j,]) > abs(meanTrue6[j]))
            pBoot5[j] = mean(abs(meanOut5[j,]) > abs(meanTrue5[j]))
            pBoot4[j] = mean(abs(meanOut4[j,]) > abs(meanTrue4[j]))
            pBoot3[j] = mean(abs(meanOut3[j,]) > abs(meanTrue3[j]))
            pBoot2[j] = mean(abs(meanOut2[j,]) > abs(meanTrue2[j]))
            pBoot1[j] = mean(abs(meanOut1[j,]) > abs(meanTrue1[j]))
            
          }
          pBootFDR6=p.adjust(pBoot6,method = "fdr",)
          pBootFDR5=p.adjust(pBoot5,method = "fdr",)
          pBootFDR4=p.adjust(pBoot4,method = "fdr",)
          pBootFDR3=p.adjust(pBoot3,method = "fdr",)
          pBootFDR2=p.adjust(pBoot2,method = "fdr",)
          pBootFDR1=p.adjust(pBoot1,method = "fdr",)
          
          
          pBoot=data.frame(pBoot1, pBoot2, pBoot3, pBoot4, pBoot5, pBoot6)
          pBootFDR=data.frame(pBootFDR1, pBootFDR2, pBootFDR3, pBootFDR4, pBootFDR5, pBootFDR6)
          
          rownames(pBoot)=rownames(geneSetEnrich)
          rownames(pBootFDR)=rownames(geneSetEnrich)
          #End permutation
          
          pBootMic=pBoot
          pBootFDRMic=pBootFDR
          
        }
      
      if(geneSetNumber=="1"){
        save_pBootMic1=data.frame(pBootMic,pBootFDRMic)
      }else{
        save_pBootMic2=data.frame(pBootMic,pBootFDRMic)
      }
      
      
      #plotting neuron genes
      heatmap3(DataGSVA.microglia, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
               col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
               Rowv=NA, Colv=NA,  scale="row",
               labRow = NA, labCol = NA, main="microglia genes") 
      
      #Compute p-values and FDR q-values for each gene set.
      # pValsTukey=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
      # pANOVA=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
      numGenes=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)
      
      for (i in 1:dim(geneSetEnrich)[1])
      { #p=t.test(geneSetEnrichSort[i,allColorGroups=="Blue"], geneSetEnrichSort[i,allColorGroups=="Purple"])
        #pVals[i]=p$p.value
        
        indGetGS=which(match(names(microGeneSet),row.names(geneSetEnrich)[i])==1)
        numGenes[i]=length(unlist(microGeneSet[indGetGS]))
        
        #ANOVA
        #   dataANOVA=data.frame(geneSetEnrichSort[i,], allColorGroups)
        #   colnames(dataANOVA)=c("GSVA_ES","color")
        #   AOVObject=aov(GSVA_ES ~ color, data = dataANOVA)
        #   resultANOVA=summary(AOVObject)
        #   pANOVA[i]=resultANOVA[[1]][["Pr(>F)"]][1]
        #   
        #   ANOVATukey=TukeyHSD(AOVObject,conf.level = 0.95)
        #   pValsTukey[i]=ANOVATukey$color[which(rownames(ANOVATukey$color)=="skyblue4-orangered4"),4] #Compare blue vs purple
      }
      
      # FDRANOVA=p.adjust(pANOVA,method = "fdr")
      # qFDRTukey=p.adjust(pValsTukey,method = "fdr")
      geneSets_Stats=data.frame(rownames(geneSetEnrichSort),pBootMic,pBootFDRMic, numGenes);
      write.xlsx(geneSets_Stats, file=paste("geneSetStats_Microglia_FXN",geneSetNumber,"_norm.xlsx",sep=''), sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)
      
      
      
      #Plot bar graphs and gene set heatmaps
      
      meanBluevsOR=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) #C2 v AD1
      meanBluevsOR4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v AD2
      meanBluevsPapaya=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v MCI1
      meanBluevsPeach=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v MCI2
      meanBluevsBlue4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) # C2 v C1
      meanPeachvsPapaya =matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1) #MCI2 v MCI1
      
      for (indGeneSetPlot in 1:dim(geneSetEnrich)[1])
      {
        dataBar=data.frame(as.data.frame(geneSetEnrich[indGeneSetPlot,indSort]), as.data.frame(allColorGroups))
        colnames(dataBar)=c("ES","colors")
        
        #Summarize stats for each group color
        statsOut = summarySE(dataBar, measurevar="ES", groupvars="colors")
        
        meanBluevsOR[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered"] # AD1 - C2
        meanBluevsOR4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered4"] #AD2 - C2
        meanBluevsPapaya[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="papayawhip"] #C2 v MC1
        meanBluevsPeach[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="peachpuff3"] #C2 vs MC2
        meanBluevsBlue4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="skyblue4"] # C2 vs C1
        meanPeachvsPapaya[indGeneSetPlot]=statsOut$ES[statsOut$color=="peachpuff3"]-statsOut$ES[statsOut$color=="papayawhip"] # MCI2 v MCI1
        
        
        
        #Re-order the colors to match the heatmaps
        statsOut$colors = factor(statsOut$colors, levels = c("skyblue4","skyblue","papayawhip","peachpuff3","orangered","orangered4"))
        
        group.colors <- c(skyblue4 = "skyblue4", skyblue = "skyblue", papayawhip="papayawhip",peachpuff3="peachpuff3", orangered ="orangered", orangered4 = "orangered4")
        
        pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia/geneSets",indGeneSetPlot,".png",sep='')
        png(pngPath, width=5,height=4,units="in",res=600)
        print({
          ggplot(statsOut, aes(x=colors, y=ES, fill=colors)) + 
            geom_bar(position=position_dodge(), stat="identity", colour="black") +
            geom_errorbar(aes(ymin=ES-se, ymax=ES+se),
                          width=.2, size=1.5,                    # Width of the error bars
                          position=position_dodge(.9))+
            theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                  text = element_text(size=28),
                  axis.text.x = element_text(face = "plain", color = "black", 
                                             size = 24),
                  axis.text.y = element_text(face = "plain", color = "black", 
                                             size = 24),
                  legend.position = "none")+
            xlab("")+
            ylab("Enrichment Score")+
            ggtitle(rownames(geneSetEnrich)[indGeneSetPlot])+
            scale_fill_manual(values=group.colors)+
            ylim(-0.5,0.5)+
            scale_x_discrete(labels=c("C1", "C2", "M1", "M2", "A1", "A2"))+
            geom_hline(yintercept=0)
          
        })
        dev.off()
        
        
        
        # #plot heatmaps for each gene set
        # 
        # #Get index of gene set associated with enrichment
        # indGetGS=which(match(names(microGeneSet),row.names(geneSetEnrich)[indGeneSetPlot])==1)
        # 
        # 
        # indGS_inData=match(as.vector(unlist(microGeneSet[indGetGS])),geneSortClean)
        # indGS_inData=indGS_inData[!is.na(indGS_inData)]
        # 
        # #FXN_Annotation_Neurons=read_excel("../Functional Categorization/20190808_NeuronFunctionalCategorization_RCluster.xlsx", sheet = "Sheet1")
        # FXN_Annotation_Neurons=AnnotationMicro.2
        # 
        # 
        # if (length(indGS_inData)>1){
        # 
        #   ## C2 v AD1 ##
        #   pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD1_rev/geneSetHeatMap_C2vAD1_",indGeneSetPlot,".png",sep='')
        #   png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #   #png(pngPath, height=700, width=650, pointsize=12)
        #   print({
        # 
        #     blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale))
        #     rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
        #     blah=blah[!is.na(blah[,1]),]
        # 
        #     indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
        # 
        #     rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
        # 
        # 
        #     hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
        #     hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
        # 
        #     #RedGreenDiff=rowMeans(blah[,indRed])-rowMeans(blah[,indGreen])
        #     RedGreenDiff=rowMeans(blah[,indGreen])-rowMeans(blah[,indRed])
        # 
        # 
        #     indBlah.RG=sort(RedGreenDiff, index.return=TRUE)$ix
        # 
        #     heatmap3(blah[indBlah.RG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
        #              col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5),
        #              Rowv=NA, Colv=NA,  scale="row",
        #              labRow = rowLabels[indBlah.RG], labCol = NA, cexRow=0.15,
        #              main=rownames(geneSetEnrich)[indGeneSetPlot])
        #   })
        #   dev.off()
        # 
        #   ## C2 v AD2 ##
        #   pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD2_rev/geneSetHeatMap_C2vAD2_",indGeneSetPlot,".png",sep='')
        #   png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #   #png(pngPath, height=700, width=650, pointsize=12)
        #   print({
        # 
        #     blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale))
        #     rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
        #     blah=blah[!is.na(blah[,1]),]
        # 
        #     indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
        # 
        #     rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
        # 
        # 
        #     hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
        #     hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
        # 
        #     PurpleGreenDiff=rowMeans(blah[,indGreen])-rowMeans(blah[,indPurple])
        # 
        #     indBlah.PG=sort(PurpleGreenDiff, index.return=TRUE)$ix
        # 
        #     heatmap3(blah[indBlah.PG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
        #              col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5),
        #              Rowv=NA, Colv=NA,  scale="row",
        #              labRow = rowLabels[indBlah.PG], labCol = NA, cexRow=0.15,
        #              main=rownames(geneSetEnrich)[indGeneSetPlot])
        #   })
        #   dev.off()
        # 
        #   ## C2 v MCI1 ##
        #   pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI1_rev/geneSetHeatMap_C2vMCI1_",indGeneSetPlot,".png",sep='')
        #   png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #   #png(pngPath, height=700, width=650, pointsize=12)
        #   print({
        # 
        #     blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale))
        #     rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
        #     blah=blah[!is.na(blah[,1]),]
        # 
        #     indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
        # 
        #     rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
        # 
        # 
        #     hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
        #     hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
        # 
        #     PapayaGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indPlum])
        # 
        #     indBlah.PapG =sort(PapayaGreenDiff, index.return=TRUE)$ix
        # 
        #     heatmap3(blah[indBlah.PapG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
        #              col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5),
        #              Rowv=NA, Colv=NA,  scale="row",
        #              labRow = rowLabels[indBlah.PapG], labCol = NA, cexRow=0.15,
        #              main=rownames(geneSetEnrich)[indGeneSetPlot])
        #   })
        #   dev.off()
        # 
        #   ## C2 v MCI2 ##
        #   pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI2_rev/geneSetHeatMap_C2vMCI2_",indGeneSetPlot,".png",sep='')
        #   png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #   #png(pngPath, height=700, width=650, pointsize=12)
        #   print({
        # 
        #     blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale))
        #     rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
        #     blah=blah[!is.na(blah[,1]),]
        # 
        #     indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
        # 
        #     rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
        # 
        # 
        #     hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
        #     hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
        # 
        #     IndianRedGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indIndianred])
        # 
        #     indBlah.IRG =sort(IndianRedGreenDiff, index.return=TRUE)$ix
        # 
        #     heatmap3(blah[indBlah.IRG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
        #              col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5),
        #              Rowv=NA, Colv=NA,  scale="row",
        #              labRow = rowLabels[indBlah.IRG], labCol = NA, cexRow=0.15,
        #              main=rownames(geneSetEnrich)[indGeneSetPlot])
        #   })
        #   dev.off()
        # 
        #   ## C2 v C1 ##
        #   pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".png",sep='')
        #   png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #   #png(pngPath, height=700, width=650, pointsize=12)
        #   print({
        # 
        #     blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale))
        #     rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
        #     blah=blah[!is.na(blah[,1]),]
        # 
        #     indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
        # 
        #     rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
        # 
        # 
        #     hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
        #     hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
        # 
        #     BlueGreenDiff= rowMeans(blah[,indGreen])-rowMeans(blah[,indBlue])
        # 
        #     indBlah.BG =sort(BlueGreenDiff, index.return=TRUE)$ix
        # 
        #     heatmap3(blah[indBlah.BG,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
        #              col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5),
        #              Rowv=NA, Colv=NA,  scale="row",
        #              labRow = rowLabels[indBlah.BG], labCol = NA, cexRow=0.15,
        #              main=rownames(geneSetEnrich)[indGeneSetPlot])
        #   })
        #   dev.off()
        # 
        #   ## MCI2 v MCI1 ##
        #   pngPath=paste("./FXN", geneSetNumber,"_norm_Neurons_MCI2vMCI1_rev/geneSetHeatMap_MCI2vMCI1_",indGeneSetPlot,".png",sep='')
        #   png(pngPath,width=3.25,height=3.25,units="in",res=600)
        #   #png(pngPath, height=700, width=650, pointsize=12)
        #   print({
        # 
        #     blah=t(apply(yDataforGSVA[indGS_inData,indSort], 1, scale))
        #     rownames(blah)=rownames(yDataforGSVA[indGS_inData,indSort])
        #     blah=blah[!is.na(blah[,1]),]
        # 
        #     indFXN2_in_GS=match(rownames(blah), FXN_Annotation_Neurons$GeneSymbol)
        # 
        #     rowLabels=paste(geneSortClean[indGS_inData],"_",FXN_Annotation_Neurons$FXN2[indFXN2_in_GS], sep = "")
        # 
        # 
        #     hrGSVA = hclust(as.dist(1-cor(t(blah), method="pearson")), method="complete")
        #     hrGSVA= hclust(dist((blah),method = "euclidean"), method="ward.D2")
        # 
        #     IRPlumDiff= rowMeans(blah[,indIndianred])-rowMeans(blah[,indPlum])
        # 
        #     indBlah.IRP =sort(IRPlumDiff, index.return=TRUE)$ix
        # 
        #     heatmap3(blah[indBlah.IRP,], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
        #              col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5),
        #              Rowv=NA, Colv=NA,  scale="row",
        #              labRow = rowLabels[indBlah.IRP], labCol = NA, cexRow=0.15,
        #              main=rownames(geneSetEnrich)[indGeneSetPlot])
        #   })
        #   dev.off()
        # 
        # 
        # 
        # }
        # graphics.off()
        # 
        # RedGreenDiffOut=cbind(RedGreenDiff[indBlah.RG[length(indBlah.RG):1]],rowLabels[indBlah.RG[length(indBlah.RG):1]])
        # PurpleGreenDiffOut=cbind(PurpleGreenDiff[indBlah.PG[length(indBlah.PG):1]],rowLabels[indBlah.PG[length(indBlah.PG):1]])
        # PapayaGreenDiffOut=cbind(PapayaGreenDiff[indBlah.PapG[length(indBlah.PapG):1]],rowLabels[indBlah.PapG[length(indBlah.PapG):1]])
        # IndianRedGreenDiffOut=cbind(IndianRedGreenDiff[indBlah.IRG[length(indBlah.IRG):1]],rowLabels[indBlah.IRG[length(indBlah.IRG):1]])
        # BlueGreenDiffOut=cbind(BlueGreenDiff[indBlah.BG[length(indBlah.BG):1]],rowLabels[indBlah.BG[length(indBlah.BG):1]])
        # IRPlumDiffOut=cbind(IRPlumDiff[indBlah.IRP[length(indBlah.IRP):1]],rowLabels[indBlah.IRP[length(indBlah.IRP):1]])
        # 
        # write.xlsx(RedGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD1_rev/geneSetHeatMap_C2vAD1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        # write.xlsx(PurpleGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vAD2_rev/geneSetHeatMap_C2vAD2_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        # write.xlsx(PapayaGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI1_rev/geneSetHeatMap_C2vMCI1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        # write.xlsx(IndianRedGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vMCI2_rev/geneSetHeatMap_C2vMCI2_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        # write.xlsx(BlueGreenDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)
        # write.xlsx(IRPlumDiffOut, file=paste("./FXN", geneSetNumber,"_norm_Neurons_C2vC1_rev/geneSetHeatMap_C2vC1_",indGeneSetPlot,".xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)

      }
      
      
      #Make ES difference bar plot for neurons
      
      
      ### C2 v AD1
      geneSetDifferenceData.C2vA1=data.frame(meanBluevsOR,rownames(geneSetEnrich),pBootFDRMic[,4])
      
      colnames(geneSetDifferenceData.C2vA1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      #geneSetDifferenceData.C2vA1[which(geneSetDifferenceData.C2vA1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia_C2vAD1_rev/ES_diff_FXN_Neu_C2vAD1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vA1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-A1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRMic[,4]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vA1)[1]))
        
      })
      dev.off()
      
      ### C2 v AD2
      geneSetDifferenceData.C2vA2=data.frame(meanBluevsOR4,rownames(geneSetEnrich),pBootFDRMic[,5])
      
      colnames(geneSetDifferenceData.C2vA2)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      #geneSetDifferenceData.C2vA2[which(geneSetDifferenceData.C2vA2$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia_C2vAD2_rev/ES_diff_FXN_Neu_C2vAD2_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vA2, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-A2)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRMic[,5]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-.5,1), xlim=c(1,dim(geneSetDifferenceData.C2vA2)[1]))
        
      })
      dev.off()
      
      ### C2 v MCI1
      geneSetDifferenceData.C2vMCI1=data.frame(meanBluevsPapaya,rownames(geneSetEnrich),pBootFDRMic[,2])
      
      colnames(geneSetDifferenceData.C2vMCI1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      #geneSetDifferenceData.C2vMCI1[which(geneSetDifferenceData.C2vMCI1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia_C2vMCI1_rev/ES_diff_FXN_Neu_C2vMCI1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vMCI1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                  fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-MCI1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRMic[,2]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vMCI1)[1]))
        
      })
      dev.off()
      
      ### C2 v MCI2
      geneSetDifferenceData.C2vMCI2=data.frame(meanBluevsPeach,rownames(geneSetEnrich),pBootFDRMic[,3])
      
      colnames(geneSetDifferenceData.C2vMCI2)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      #geneSetDifferenceData.C2vMCI2[which(geneSetDifferenceData.C2vMCI2$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia_C2vMCI2_rev/ES_diff_FXN_Neu_C2vMCI2_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vMCI2, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                  fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          #scale_fill_manual(name = "area", values=c("red","grey50"))+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-MCI2)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRMic[,3]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vMCI2)[1]))
        
      })
      dev.off()
      
      ### C2 v C1
      geneSetDifferenceData.C2vC1=data.frame(meanBluevsBlue4,rownames(geneSetEnrich),pBootFDRMic[,1])
      
      colnames(geneSetDifferenceData.C2vC1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      #geneSetDifferenceData.C2vC1[which(geneSetDifferenceData.C2vC1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia_C2vC1_rev/ES_diff_FXN_Neu_C2vC1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.C2vC1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (C2-C1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRMic[,1]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]))
        
      })
      dev.off()
      
      ### MCI2 v MCI1
      geneSetDifferenceData.MCI2vMCI1=data.frame(meanPeachvsPapaya,rownames(geneSetEnrich),pBootFDRMic[,6])
      
      colnames(geneSetDifferenceData.MCI2vMCI1)=c("difference","geneSetName", 'pBootFDR')
      
      #Fix the spelling of "miscellaneous
      geneSetDifferenceData.MCI2vMCI1[which(geneSetDifferenceData.MCI2vMCI1$geneSetName=="Miscelaneous"),2]="Miscellaneous"
      
      
      
      pngPath=paste("./FXN",geneSetNumber,"_norm_Microglia_MCI2vMCI1_rev/ES_diff_FXN_Neu_MCI2vMCI1_",geneSetNumber,".pdf",sep='')
      #png(pngPath, width=6,height=4,units="in",res=600)
      pdf(pngPath, width=6,height=4)
      
      print({
        ggplot(geneSetDifferenceData.MCI2vMCI1, aes(x=reorder(geneSetName,-difference), y=difference, 
                                                    fill=factor(ifelse(pBootFDR<=0.05,"Highlighted","Normal")))) + 
          
          geom_bar(position=position_dodge(), stat="identity", colour="black") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
                text = element_text(size=20),
                axis.text.x = element_text(face = "plain", color = "black", 
                                           size = 18),
                axis.text.y = element_text(face = "plain", color = "black", 
                                           size = 16),
                legend.position = "none")+
          scale_fill_manual(name = "area", values=c("red","grey50"))+
          xlab("Gene Set")+
          ylab("ES Diff. (MCI2-MCI1)")+
          geom_hline(yintercept=0)+
          geom_text(aes(label=format((pBootFDRMic[,6]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
          geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
          coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]))
        
      })
      dev.off()   
      
      
      
      
    
    
    
    
    
    
    
    

    
    
    
    #----Make color bars----
    color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
      scale = (length(lut)-1)/(max-min)
      
      #dev.new(width=1.75, height=5)
      plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
      axis(2, ticks, las=1)
      for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
      }
    }
    
    
    #----Make color bars----
    color.bar.log <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
      scale = (length(lut)-1)/(max-min)
      
      #dev.new(width=1.75, height=5)
      plot(c(0,1), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title, log='y')
      axis(2, ticks, las=1)
      for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
      }
    }
    
    

    
    

# ---- Write the data out ----
dataClustOut=data.frame(geneLabelsClust,clusterNumber, tau, cellTypeMax,yMatClust);
colnames(dataClustOut)=c("GeneSymbol", "ClusterNumber","Tau", "TauCellType",colnames(y));
write.xlsx(dataClustOut, file="ClusteredGeneTypeData.xlsx", sheetName="ClusteredGeneTypeData", col.names=TRUE, row.names=FALSE, append=FALSE)

