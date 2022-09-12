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
#library("DESeq2")
library("ggplot2")
#library("scater")
library("ggpubr")
library("GSVA")
library("RColorBrewer")
library("vsn")
library("limma")
library("Hmisc")
library("plotrix") 
library(clusterProfiler)
library(AnnotationDbi)

source("../Custom R Functions/Mode.R")
source("../Custom R Functions/summarySE.R")
source("../Custom R Functions/rcorr.adjustFDR.R")

library(factoextra)
library(cluster)
library("enrichplot")
library("rrvgo")

#Set heatmap3 color bar parameters
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red")

breakBarColorsCor=c(-200,seq(-1, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
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

png("../CTDCluster.png", width=7,height=7,units="in",res=600)
par(mar=c(0,0,0,0))
heatmap3(yMatCTDZ,  
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=as.dendrogram(hrCTD), Colv=as.dendrogram(hcCTD),  scale="none",
         RowSideColors=mycolhcCTD, RowSideLabs=NA, cexRow=36,cexCol =1,
         labRow = NA, labCol = colnames(yCTD), 
         margins = c(5,25), 
         ColSideWidth=.4)
dev.off()

rowOrder=order.dendrogram(t1CTD);
rowOrderCol=order.dendrogram(as.dendrogram(hcCTD));


#Prepare clustered data for writing to table
yMatClust=yMatCTD[rowOrder,rowOrderCol];
geneLabelsClustCTD=geneLabelsCTD[rowOrder];
clusterNumberCTD=myclCTD[rowOrder];


#write.xlsx(dataClustOut, file="xxxlsx", sheetName="xx", col.names=TRUE, row.names=FALSE, append=FALSE)



# ---- Read and Plot Mt Sinai Data ----

mtSinaiData=read.csv("./Data Pre-Processing/Mayo_Final.csv", sep=",")

y=mtSinaiData[2:dim(mtSinaiData)[1],4:dim(mtSinaiData)[2]]; #Careful of indexing
geneLabelsMtSinai=mtSinaiData$geneSymbSort[2:nrow(mtSinaiData)];    #Careful of indexing 
phenoMtSinai=(mtSinaiData[1,4:dim(mtSinaiData)[2]])

phenoMtSinai <- data.frame(lapply(phenoMtSinai, as.character), stringsAsFactors=FALSE)
phenoMtSinai[phenoMtSinai=="Control"]="CTRL"

yMat=apply(matrix(unlist(y),ncol=dim(y)[2],byrow=FALSE),2,as.numeric)


colorBar=ifelse(phenoMtSinai=="CTRL","navy","red3")

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
ySortZAD=ySortZ[,phenoMtSinai=="AD"];

hcSortCTRL= hclust(dist(t(ySortZCTRL),method = "euclidean"), method="ward.D2")
heatmap3(ySortZCTRL, ColSideColors =c(colorBar)[1:sum(phenoMtSinai=="CTRL")], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortCTRL),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 

hcSortAD= hclust(dist(t(ySortZAD),method = "euclidean"), method="ward.D2")
heatmap3(ySortZAD, ColSideColors =c(colorBar)[sum(phenoMtSinai=="CTRL")+1:length(phenoMtSinai)], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=as.dendrogram(hcSortAD),  scale="none",
         labRow = NA, labCol = NA, 
         margins = c(5,20)) 


#----Cell Type Normalization----
#Setup phenotype used for DEseq2 format below
phenoNum=as.integer(ifelse(phenoMtSinai=="CTRL",0,1))
condition=as.factor(c(phenoNum))

#clean cases 
indGenesKeep <- complete.cases(ySort)
ySortClean=ySort[indGenesKeep,]
geneSortClean=geneSort[indGenesKeep]
clusterNumberCTDClean=clusterNumberCTD[indGenesKeep]  

# #Clean up data - keep rows that are not NA
# indGenesKeep=which(!is.na(ySort[,1]))
# ySortClean=ySort[indGenesKeep,]
# geneSortClean=geneSort[indGenesKeep]
# clusterNumberCTDClean=clusterNumberCTD[indGenesKeep]  


#Initialize normalized data
ySortNorm=ySortClean

#create matrix with gene annotation for later extraction of genes in various subgroups
ySortNorm_ann <- ySortClean 
rownames(ySortNorm_ann) <- geneSortClean

#Normalize neurons - SVN approach
ySortNormNeuronGenes=normalizeVSN(ySortNorm[clusterNumberCTDClean==13,])
ySortNeuronGenes=ySortClean[clusterNumberCTDClean==13,]
ySortNorm[clusterNumberCTDClean==13,]=normalizeVSN(ySortNorm[clusterNumberCTDClean==13,])

neuroGeneNames <- rownames(ySortNorm_ann[clusterNumberCTDClean==13,])

#Normalize astrocytes - SVN approach
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


#AD
hcADCut = cutree(hcSortAD, k=2); 
colorhcAD = rainbow(length(unique(hcADCut)), start=0.3, end=0.6); 
colorhcAD=c("orangered","orangered4")
colorhcAD = colorhcAD[as.vector(hcADCut)]

startColorIndAD=sum(phenoMtSinai=="CTRL")+1
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
combinedColorALL=cbind(c(colorBar[startColorIndCTRL:endColorIndCTRL], colorBar[startColorIndAD:endColorIndAD]), 
                       c(colorhcCTRL,colorhcAD))

allColorGroups=c(colorhcCTRL,colorhcAD)

indBlue=which(c(colorhcCTRL,colorhcAD)=="skyblue4")
indGreen=which(c(colorhcCTRL,colorhcAD)=="skyblue")
indRed=which(c(colorhcCTRL,colorhcAD)=="orangered")
indPurple=which(c(colorhcCTRL,colorhcAD)=="orangered4")

indSort=c(indBlue,indGreen,indRed,indPurple)
allColorGroups=allColorGroups[indSort]
#geneSetEnrichSort=geneSetEnrich[,indSort]

#Write out neuron genes sorted with samples sorted by indSort
yNeuMtSinai=ySortNormNeuronGenes[,indSort]
rownames(yNeuMtSinai)=geneSortClean[clusterNumberCTDClean==13]

#Write out astrocyte genes sorted with samples sorted by indSort
yAstMayo=ySortNormAstrocyteGenes[,indSort]
rownames(yAstMayo)=geneSortClean[clusterNumberCTDClean==4]

#save(yAstMayo, phenoMayo, file="../ADGeneSignature/SigMayo.RData")

#Write out microglia genes sorted with samples sorted by indSort
yMicMtSinai=ySortNormMicrogliaGenes[,indSort]
rownames(yMicMtSinai)=geneSortClean[clusterNumberCTDClean==17]

phenoMayo=allColorGroups


#----FC----
#Compute histogram comparing fold change between normlized and non-normlized genes

#Picking groups 1 for CTRL and AD
indCTRL_FC=which(phenoMtSinai=="CTRL" & c(hcCTRLCut==1,rep(NA,endColorIndAD-startColorIndAD+1)))
indAD_FC=which(phenoMtSinai=="AD" & c(rep(NA,endColorIndCTRL), hcADCut==1))


#Compute histogram based on normalized genes
#Neurons


ySortNeuronGenesZ=t(apply(ySortNeuronGenes,1,scale)); #z-score the data
ySortNormNeuronGenesZ=t(apply(ySortNormNeuronGenes,1,scale)); #z-score the data


meanCTRL=apply(ySortNeuronGenesZ[,indBlue],1,mean)
meanAD=apply(ySortNeuronGenesZ[,indPurple],1,mean)
logFC1=meanAD-meanCTRL


meanCTRL=apply(ySortNormNeuronGenesZ[,indCTRL_FC],1,mean)
meanAD=apply(ySortNormNeuronGenesZ[,indAD_FC],1,mean)
logFC2=meanAD-meanCTRL

neuoronLogFC=data.frame(
  logFC=c(logFC1,logFC2),
  norm=rep(c("Not Normalized","Normalized"),each=length(logFC1)))

pngPath=paste("Normalization__neurons_hist.png",sep='')
png(pngPath, width=4,height=4,units="in",res=600, pointsize = 20)
ggdensity(neuoronLogFC, x="logFC", y = "..density..",
          color="norm", fill = "norm", palette = "jco",
          size=1, )+labs(title="Neuron Genes")+xlab("Mean(A2-C1)")+ylab("Density")+xlim(-3,3)+
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
logFC1=meanAD-meanCTRL


meanCTRL=apply(ySortNormAstrocyteGenesZ[,indCTRL_FC],1,mean)
meanAD=apply(ySortNormAstrocyteGenesZ[,indAD_FC],1,mean)
logFC2=meanAD-meanCTRL


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

microgliaLogFC=data.frame(
  logFC=c(logFC1,logFC2),
  norm=rep(c("Not Normalized","Normalized"),each=length(logFC1)))


png("Normalization_microglia_hist.png", width=4,height=4,units="in",res=600, pointsize = 20)

ggdensity(microgliaLogFC, x="logFC", y = "..density..",
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





#----Normalization Plotting---- 
#Put everything in same heatmap
ySortNormZ=t(apply(ySortNorm,1,scale)); #z-score the data
ySortZNeurons=t(apply(ySortNormNeuronGenes,1,scale)); #z-score the data
ySortZAstrocytes=t(apply(ySortNormAstrocyteGenes,1,scale)); #z-score the data
ySortZMicroglia=t(apply(ySortNormMicrogliaGenes,1,scale)); #z-score the data


# ---- Unsupervised self-clustering of Cell Type specific gene groups --- #

###NEURONS

clust.neurons_avg <- hclust(dist((ySortZNeurons), method = "euclidean"), method ="average")
clust.neurons_ward <- hclust(dist((ySortZNeurons), method = "euclidean"), method ="ward.D2")
cort.neu <- cor(t(ySortZNeurons), method="spearman")
clust.neurons <- hclust(as.dist(1-cort.neu),method = "complete")
#clust.neurons <- eclust(ySortZNeurons, "hclust", k=NULL, hc_metric = "spearman", hc_method = "complete", graph = FALSE )

#plot(as.dendrogram(clust.neurons))

# gap statistics for internal cluster number validation
# png("Gap statistics neurons.png", width=4,height=4,units="in",res=600, pointsize = 20)
# 
# fviz_nbclust(ySortZNeurons, hcut, method = "gap_stat", diss = as.dist(1-cort.neu), k.max = 5)
# 
# dev.off()

#cut dendrogram into indicated clusters
clust.neurons.cut <- cutree(clust.neurons, k=5)

colorBrewer.clust.neu =colorRampPalette(brewer.pal(8,"Spectral"))
mycolclust.neu =  sample(colorBrewer.clust.neu(length(unique(clust.neurons.cut))))
mycolclust.neu = mycolclust.neu[as.vector(clust.neurons.cut)]

#png("Norm Neuron Genes self clustering (correlation) k5.png", width=4,height=4,units="in",res=600, pointsize = 20)

heatmap3(ySortZNeurons[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv= as.dendrogram(clust.neurons), Colv=NA,  scale="none",
         RowSideColors = mycolclust.neu, RowSideLabs=NA,
         labRow = NA, labCol = NA, main="Normalized Neuron Genes (celltype individual clustering)")

#dev.off()

# A - clust 2 (339) - M1
# B - clust 4 (167) - M2
# C - clust 3 (83) - M3
# D - clust 5 (246) - M4
# E - clust 1 (188) - M5



### ASTROCYTES ###


clust.astrocytes_avg <- hclust(dist((ySortZAstrocytes), method = "euclidean"), method ="average")
clust.astrocytes_ward <- hclust(dist((ySortZAstrocytes), method = "euclidean"), method ="ward.D2")
cort.ast <- cor(t(ySortZAstrocytes), method = "spearman")
clust.astrocytes <- hclust(as.dist(1-cort.ast),method = "complete")
#clust.astrocytes <- eclust(ySortZAstrocytes, "hclust", k=NULL, hc_metric = "spearman", hc_method = "complete", graph = FALSE )

# #gap statistics for cluster number validation
# png("Gap statistics Astrocytes.png", width=4,height=4,units="in",res=600, pointsize = 20)
# 
# fviz_nbclust(ySortZAstrocytes, hcut, method = "gap_stat", diss = as.dist(1-cort.ast), k.max = 5)
# 
# dev.off()

#plot(as.dendrogram(clust.astrocytes))

clust.astrocytes.cut <- cutree(clust.astrocytes, k = 5)

colorBrewer.clust.ast =colorRampPalette(brewer.pal(8,"Spectral"))
mycolclust.ast =  sample(colorBrewer.clust.ast(length(unique(clust.astrocytes.cut))))
mycolclust.ast = mycolclust.ast[as.vector(clust.astrocytes.cut)]

#png("Norm Astrocyte Genes self clustering (correlation) k5.png", width=4,height=4,units="in",res=600, pointsize = 20)

heatmap3(ySortZAstrocytes[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv= as.dendrogram(clust.astrocytes), Colv=NA,  scale="none",
         RowSideColors = mycolclust.ast, RowSideLabs=NA,
         labRow = NA, labCol = NA, main="Normalized Astrocyte Genes (celltype individual clustering)")

#dev.off()

# A - clust 4 (330) - M1 
# B - clust 3 (174) - M2
# C - clust 5 (249) - M3
# D - clust 2 (361) - M4
# E - clust 1 (355) - M5

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

# clust.microglia_avg <- hclust(dist((ySortZMicroglia), method = "euclidean"), method ="average")
# clust.microglia_ward <- hclust(dist((ySortZMicroglia), method = "euclidean"), method ="ward.D2")
cort.mic <- cor(t(ySortZMicroglia), method="spearman")
clust.microglia <- hclust(as.dist(1-cort.mic),method = "complete")

# clust.microglia <- eclust(ySortZMicroglia, "hclust", k=NULL, hc_metric = "spearman", hc_method = "complete", graph = FALSE )

#gap statistics for cluster number validation
# png("Gap statistics microglia.png", width=4,height=4,units="in",res=600, pointsize = 20)
# 
# fviz_nbclust(ySortZMicroglia, hcut, method = "gap_stat", diss = as.dist(1-cort), k.max = 5)
# 
# dev.off()

#plot(as.dendrogram(clust.microglia))

clust.microglia.cut <- cutree(clust.microglia, k=5) 

table(clust.microglia.cut)

colorBrewer.clust.mic =colorRampPalette(brewer.pal(8,"Spectral"))
mycolclust.mic =  sample(colorBrewer.clust.mic(length(unique(clust.microglia.cut))))
mycolclust.mic = mycolclust.mic[as.vector(clust.microglia.cut)]

#png("Norm Microglia Genes self clustering (correlation) k5.png", width=4,height=4,units="in",res=600, pointsize = 20)

heatmap3(ySortZMicroglia[,indSort], ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv= as.dendrogram(clust.microglia), Colv=NA,  scale="none",
         RowSideColors = mycolclust.mic, RowSideLabs=NA,
         labRow = NA, labCol = NA, main="Normalized Microglia Genes (celltype individual clustering)")

#dev.off()

# A - clust 3 (405) - M1 
# B - clust 4 (257) - M2
# C - clust 1 (260) - M3
# D - clust 2 (277) - M4
# E - clust 5 (196) - M5

# ---- GSVA on Cell Type Normalized Data ----

#Choose the data to run GSVA on.
# yDataforGSVA=(ySortNorm);
# rownames(yDataforGSVA)=geneSortClean


# #----First, Wilcoxon test Astrocytes----
# wilcoxTest=matrix(NA,nrow = nrow(yDataforGSVA),ncol=6)
# colnames(wilcoxTest)=c("C1","SEM_C1","AD2","SEM_AD2","p","pFDR")
# rownames(wilcoxTest)=rownames(yDataforGSVA)
# for (i in 1:nrow(yDataforGSVA)) {
#   wilcoxTest[i,1]=mean(yDataforGSVA[i,indBlue])
#   wilcoxTest[i,2]=std.error(yDataforGSVA[i,indBlue])
#   wilcoxTest[i,3]=mean(yDataforGSVA[i,indPurple])
#   wilcoxTest[i,4]=std.error(yDataforGSVA[i,indPurple])
#   wilcoxTest[i,5]=wilcox.test(yDataforGSVA[i,indPurple], yDataforGSVA[i,indBlue], alternative = "two.sided")$p.value
#   
# }
# 
# wilcoxTestAstrocytes=wilcoxTest[clusterNumberCTDClean==4,]
# wilcoxTestAstrocytes[,6]=p.adjust(wilcoxTestAstrocytes[,5], method="fdr")
# write.xlsx(wilcoxTestAstrocytes, file="astrocyteWilcoxTest_Mayo.xlsx", sheetName="sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)



#Choose the gene set number
geneSetNumber="1"

#Choose whether or not to compute p values
compPermute=1




#@@@@@@@ UNANNOTATED @@@@@@# 


#----____NEURONS----

ySortZNeurons_ann <- ySortZNeurons
rownames(ySortZNeurons_ann) <- neuroGeneNames

#seperate data based on clusters
neuroDataClust.1 <- ySortZNeurons_ann[clust.neurons.cut == 1,]
neuroDataClust.2 <- ySortZNeurons_ann[clust.neurons.cut == 2,]
neuroDataClust.3 <- ySortZNeurons_ann[clust.neurons.cut == 3,]
neuroDataClust.4 <- ySortZNeurons_ann[clust.neurons.cut == 4,]
neuroDataClust.5 <- ySortZNeurons_ann[clust.neurons.cut == 5,]

#export micrglia data with cluster annotation
writeNeuroData <- as.data.frame(ySortZNeurons_ann)
writeNeuroData$symbol <- rownames(writeNeuroData)
writeNeuroData$cluster <- "na"

writeNeuroData[clust.neurons.cut == 1,]$cluster <- "cluster1"
writeNeuroData[clust.neurons.cut == 2,]$cluster <- "cluster2"
writeNeuroData[clust.neurons.cut == 3,]$cluster <- "cluster3"
writeNeuroData[clust.neurons.cut == 4,]$cluster <- "cluster4"
writeNeuroData[clust.neurons.cut == 5,]$cluster <- "cluster5"

writeNeuroData <- (writeNeuroData[,c("symbol","cluster")])

write.xlsx(writeNeuroData, file=paste("./unannotated/AstrocyteClusterGeneSet.xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)

# A - clust 2 (339) - M1
# B - clust 4 (167) - M2
# C - clust 3 (83) - M3
# D - clust 5 (246) - M4
# E - clust 1 (188) - M5

#save gene names of clusters
neuroGenesClust.1 <- rownames(ySortZNeurons_ann[clust.neurons.cut == 1,])
neuroGenesClust.2 <- rownames(ySortZNeurons_ann[clust.neurons.cut == 2,])
neuroGenesClust.3 <- rownames(ySortZNeurons_ann[clust.neurons.cut == 3,])
neuroGenesClust.4 <- rownames(ySortZNeurons_ann[clust.neurons.cut == 4,])
neuroGenesClust.5 <- rownames(ySortZNeurons_ann[clust.neurons.cut == 5,])

### GO analysis to identify Biological Functions of the clusters ###

#convert gene symbols to EntrezIDs 
entrezClust1 <- mapIds(org.Hs.eg.db, neuroGenesClust.1, "ENTREZID", "SYMBOL")
entrezClust2 <- mapIds(org.Hs.eg.db, neuroGenesClust.2, "ENTREZID", "SYMBOL")
entrezClust3 <- mapIds(org.Hs.eg.db, neuroGenesClust.3, "ENTREZID", "SYMBOL")
entrezClust4 <- mapIds(org.Hs.eg.db, neuroGenesClust.4, "ENTREZID", "SYMBOL")
entrezClust5 <- mapIds(org.Hs.eg.db, neuroGenesClust.5, "ENTREZID", "SYMBOL")

#perfrom GO 

GOclust1 <- enrichGO(entrezClust1, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # N5
GOclust2 <- enrichGO(entrezClust2, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # N1
GOclust3 <- enrichGO(entrezClust3, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # N3
GOclust4 <- enrichGO(entrezClust4, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # N2
GOclust5 <- enrichGO(entrezClust5, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # N4

#plot GO results 
png("./GO/GO_neuron_N5.png")
barplot(GOclust1)
dev.off()

png("./GO/GO_neuron_N1")
barplot(GOclust2)
dev.off()

png("./GO/GO_neuron_N3")
barplot(GOclust3)
dev.off()

png("./GO/GO_neuron_N2")
barplot(GOclust4)
dev.off()

png("./GO/GO_neuron_N4")
barplot(GOclust5)
dev.off()


#generate geneSet for GSVA
neuroGenes.total <- list(neuroGenesClust.1,neuroGenesClust.2,neuroGenesClust.3,neuroGenesClust.4,neuroGenesClust.5)
list.names <- c("N-5","N-1","N-3","N-2","N-4")
names(neuroGenes.total) <- list.names
neuroGeneSet =lapply(neuroGenes.total, function(x) x[!is.na(x)])

# DATA for GSVA in Microglia 
#DataGSVA.microglia <- yMicrogliaClean #are the normalized data (after vsn) but only of microglia (clust = 17) // NOT z-scored
DataGSVA.neuron <- ySortNormNeuronGenes #are the normalized data (after vsn) but only of microglia (clust = 17) // NOT z-scored
rownames(DataGSVA.neuron)  <- neuroGeneNames

geneSetEnrich=gsva(DataGSVA.neuron, neuroGeneSet, mx.diff=TRUE, kcdf="Gaussian")

geneSetEnrich_Neuron=geneSetEnrich


geneSetEnrichSort=geneSetEnrich[,indSort] #sort samples by phenotype
geneSetEnrichSortZ=t(apply(geneSetEnrichSort,1,scale)); #z-score the data

hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Neurons_unannotated/GSVA_cluster_v2.png",sep='')
png(pngPath, width=5,height=4,units="in",res=1200, pointsize = 16)

heatmap3(geneSetEnrichSortZ, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=4.5), 
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = colnames(neuroGeneSet), labCol = NA, main="GSVA for neuron gene sets") 

dev.off()


corGSVA=cor(t(geneSetEnrichSortZ))
hCor=hclust(as.dist((1-corGSVA)/2))

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Neurons_unannotated/GSVA_corr_v2.png",sep='')
png(pngPath, width=6,height=4,units="in",res=1200, pointsize = 10)

heatmap3(corGSVA,
         col=barColorsCor, breaks=breakBarColorsCor,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",
         labRow = colnames(neuroGeneSet), labCol = NA, main="GSVA for neuron gene sets") 

dev.off()

blah=DataGSVA.neuron


#Permutation p-value
if (compPermute==1){
  R=1000
  meanOut3=matrix(NA,nrow=length(neuroGeneSet), ncol = R)
  meanOut2=matrix(NA,nrow=length(neuroGeneSet), ncol = R)
  meanOut1=matrix(NA,nrow=length(neuroGeneSet), ncol = R)
  
  for (i in 1:R)
  {
    blah2=blah
    rownames(blah2)=rownames(blah)[sample(1:dim(blah)[1])]
    gsvaBoot=gsva(blah2,neuroGeneSet, mx.diff=TRUE, kcdf="Gaussian")
    
    
    meanOut3[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPurple]) #mean C2 - Ad2
    meanOut2[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indRed])    #mean C2 - AD1
    meanOut1[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indBlue])   #mean C2 - C1
    
    
    print(i)
  }
  meanTrue3=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPurple]) #mean C2 - AD2
  meanTrue2=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indRed])   #mean C2 - AD1
  meanTrue1=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indBlue]) #mean C2 - C1
  
  pBoot3=matrix(NA,nrow=length(neuroGeneSet), ncol = 1)
  pBoot2=matrix(NA,nrow=length(neuroGeneSet), ncol = 1)
  pBoot1=matrix(NA,nrow=length(neuroGeneSet), ncol = 1)
  
  for (j in 1:length(neuroGeneSet))
  {
    pBoot3[j] = mean(abs(meanOut3[j,]) > abs(meanTrue3[j])) # AD2 - C2
    pBoot2[j] = mean(abs(meanOut2[j,]) > abs(meanTrue2[j])) # AD1 - C2
    pBoot1[j] = mean(abs(meanOut1[j,]) > abs(meanTrue1[j])) # C1 - C2
    
  }
  pBootFDR3=p.adjust(pBoot3,method = "fdr",)
  pBootFDR2=p.adjust(pBoot2,method = "fdr",)
  pBootFDR1=p.adjust(pBoot1,method = "fdr",)
  
  
  pBoot=data.frame(pBoot1, pBoot2, pBoot3)
  pBootFDR=data.frame(pBootFDR1, pBootFDR2, pBootFDR3)
  
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
heatmap3(DataGSVA.neuron, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA, main="neuron genes") 

#Compute p-values and FDR q-values for each gene set.
numGenes=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)

for (i in 1:dim(geneSetEnrich)[1])
{
  
  indGetGS=which(match(names(neuroGeneSet),row.names(geneSetEnrich)[i])==1)
  numGenes[i]=length(unlist(neuroGeneSet[indGetGS]))
  
}

geneSets_Stats=data.frame(rownames(geneSetEnrichSort),pBootNeu,pBootFDRNeu, numGenes);
write.xlsx(geneSets_Stats, file=paste("geneSetStats_Neurons_FXN_",geneSetNumber,"_norm_unannotated.xlsx",sep=''), sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)



#Plot bar graphs and gene set heatmaps

meanBluevsOR=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)
meanBluevsOR4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)
meanBluevsBlue4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)

for (indGeneSetPlot in 1:dim(geneSetEnrich)[1])
  
{
  dataBar=data.frame(as.data.frame(geneSetEnrich[indGeneSetPlot,indSort]), as.data.frame(allColorGroups))
  colnames(dataBar)=c("ES","colors")
  
  #Summarize stats for each group color
  statsOut = summarySE(dataBar, measurevar="ES", groupvars="colors") #count/mean/SD/CI 
  
  meanBluevsOR[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered"] # C2 - Ad1
  meanBluevsOR4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered4"] # C2 - AD2
  meanBluevsBlue4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="skyblue4"] # C2 - C1
  
  #Re-order the colors to match the heatmaps
  statsOut$colors = factor(statsOut$colors, levels = c("skyblue4","skyblue","orangered","orangered4"))
  
  group.colors <- c(skyblue4 = "skyblue4", skyblue = "skyblue", orangered ="orangered", orangered4 = "orangered4")
  
  pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Neurons_unannotated/geneSets",indGeneSetPlot,".png",sep='') 
  png(pngPath, width=4,height=4,units="in",res=600)
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
            axis.title.x = element_text(face = "plain", color = "black", 
                                        size = 24),
            axis.text.y = element_text(face = "plain", color = "black", 
                                       size = 24),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5))+
      xlab("")+
      ylab("Enrichment Score")+
      ggtitle(rownames(geneSetEnrich)[indGeneSetPlot])+
      scale_fill_manual(values=group.colors)+
      ylim(-1,1)+
      scale_x_discrete(labels=c("C1", "C2", "A1", "A2"))+
      geom_hline(yintercept=0)
    
  })
  dev.off()
  
}
#Make ES difference bar plot for neurons

### C2 v AD1
geneSetDifferenceData.C2vA1=data.frame(meanBluevsOR,rownames(geneSetEnrich),pBootFDRNeu[,2])
colnames(geneSetDifferenceData.C2vA1)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Neurons_C2vAD1_rev/ES_diff_FXN_Neu_C2vAD1_",geneSetNumber,".pdf",sep='')
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
    geom_text(aes(label=format((pBootFDRNeu[,2]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
    geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vA1)[1]))
  
})
dev.off()

### C2 v AD2
geneSetDifferenceData.C2vA2=data.frame(meanBluevsOR4,rownames(geneSetEnrich),pBootFDRNeu[,3])
colnames(geneSetDifferenceData.C2vA2)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Neurons_C2vAD2_rev/ES_diff_FXN_Neu_C2vAD2_",geneSetNumber,".pdf",sep='')
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
    #scale_fill_manual(name = "area", values=c("red","grey50"))+
    scale_fill_manual(name = "area", values=c("red", "grey50"))+
    xlab("Gene Set")+
    ylab("ES Diff. (C2-A2)")+
    geom_hline(yintercept=0)+
    geom_text(aes(label=format((pBootFDRNeu[,3]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
    geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vA2)[1]))
  
})
dev.off()

### C2 v C1
geneSetDifferenceData.C2vC1=data.frame(meanBluevsBlue4,rownames(geneSetEnrich),pBootFDRNeu[,1])
colnames(geneSetDifferenceData.C2vC1)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Neurons_C2vC1_rev/ES_diff_FXN_Neu_C2vC1_",geneSetNumber,".pdf",sep='')
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
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]))
  
})
dev.off()



#----____ASTROCYTES----

ySortZAstrocytes_ann <- ySortZAstrocytes
rownames(ySortZAstrocytes_ann) <- astroGeneNames

#seperate data based on clusters
astroDataClust.1 <- ySortZAstrocytes_ann[clust.astrocytes.cut == 1,]
astroDataClust.2 <- ySortZAstrocytes_ann[clust.astrocytes.cut == 2,]
astroDataClust.3 <- ySortZAstrocytes_ann[clust.astrocytes.cut == 3,]
astroDataClust.4 <- ySortZAstrocytes_ann[clust.astrocytes.cut == 4,]
astroDataClust.5 <- ySortZAstrocytes_ann[clust.astrocytes.cut == 5,]

#export micrglia data with cluster annotation
writeAstroData <- as.data.frame(ySortZAstrocytes_ann)
writeAstroData$symbol <- rownames(writeAstroData)
writeAstroData$cluster <- "na"

writeAstroData[clust.astrocytes.cut == 1,]$cluster <- "cluster1"
writeAstroData[clust.astrocytes.cut == 2,]$cluster <- "cluster2"
writeAstroData[clust.astrocytes.cut == 3,]$cluster <- "cluster3"
writeAstroData[clust.astrocytes.cut == 4,]$cluster <- "cluster4"
writeAstroData[clust.astrocytes.cut == 5,]$cluster <- "cluster5"

writeAstroData <- (writeAstroData[,c("symbol","cluster")])

write.xlsx(writeAstroData, file=paste("./unannotated/AstroClusterGeneSet.xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)

# A - clust 4 (330) - M1 
# B - clust 3 (174) - M2
# C - clust 5 (249) - M3
# D - clust 2 (361) - M4
# E - clust 1 (355) - M5

#save gene names of clusters
astroGenesClust.1 <- rownames(ySortZAstrocytes_ann[clust.astrocytes.cut == 1,])
astroGenesClust.2 <- rownames(ySortZAstrocytes_ann[clust.astrocytes.cut == 2,])
astroGenesClust.3 <- rownames(ySortZAstrocytes_ann[clust.astrocytes.cut == 3,])
astroGenesClust.4 <- rownames(ySortZAstrocytes_ann[clust.astrocytes.cut == 4,])
astroGenesClust.5 <- rownames(ySortZAstrocytes_ann[clust.astrocytes.cut == 5,])

### GO analysis to identify Biological Functions of the clusters ###

#convert gene symbols to EntrezIDs 
entrezClust1 <- mapIds(org.Hs.eg.db, astroGenesClust.1, "ENTREZID", "SYMBOL")
entrezClust2 <- mapIds(org.Hs.eg.db, astroGenesClust.2, "ENTREZID", "SYMBOL")
entrezClust3 <- mapIds(org.Hs.eg.db, astroGenesClust.3, "ENTREZID", "SYMBOL")
entrezClust4 <- mapIds(org.Hs.eg.db, astroGenesClust.4, "ENTREZID", "SYMBOL")
entrezClust5 <- mapIds(org.Hs.eg.db, astroGenesClust.5, "ENTREZID", "SYMBOL")

#perfrom GO 

GOclust1 <- enrichGO(entrezClust1, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # A-4
GOclust2 <- enrichGO(entrezClust2, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # A-5
GOclust3 <- enrichGO(entrezClust3, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # A-3
GOclust4 <- enrichGO(entrezClust4, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # A-1
GOclust5 <- enrichGO(entrezClust5, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.05 , pAdjustMethod = "fdr", readable= TRUE) # A-2

#plot GO results 
png("./GO/GO_astro_A5.png")
barplot(GOclust1)
dev.off()

png("./GO/GO_astro_A4.png")
barplot(GOclust2)
dev.off()

png("./GO/GO_astro_A2.png")
barplot(GOclust3)
dev.off()

png("./GO/GO_astro_A1.png")
barplot(GOclust4)
dev.off()

png("./GO/GO_astro_A3.png")
barplot(GOclust5)
dev.off()

#generate geneSet for GSVA
astroGenes.total <- list(astroGenesClust.1,astroGenesClust.2,astroGenesClust.3,astroGenesClust.4,astroGenesClust.5)
list.names <- c("A-5","A-4","A-2","A-1","A-3")
names(astroGenes.total) <- list.names
astroGeneSet =lapply(astroGenes.total, function(x) x[!is.na(x)])

# DATA for GSVA in Microglia 
#DataGSVA.microglia <- yMicrogliaClean #are the normalized data (after vsn) but only of microglia (clust = 17) // NOT z-scored
DataGSVA.astrocytes <- ySortNormAstrocyteGenes #are the normalized data (after vsn) but only of microglia (clust = 17) // NOT z-scored
rownames(DataGSVA.astrocytes)  <- astroGeneNames

geneSetEnrich=gsva(DataGSVA.astrocytes, astroGeneSet, mx.diff=TRUE, kcdf="Gaussian")

geneSetEnrich_Astrocyte=geneSetEnrich


geneSetEnrichSort=geneSetEnrich[,indSort] #sort samples by phenotype
geneSetEnrichSortZ=t(apply(geneSetEnrichSort,1,scale)); #z-score the data

hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_AStrocytes_unannotated/GSVA_cluster_v2.png",sep='')
png(pngPath, width=5,height=4,units="in",res=1200, pointsize = 16)

heatmap3(geneSetEnrichSortZ, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=4.5), 
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = colnames(astroGeneSet), labCol = NA, main="GSVA for neuron gene sets") 

dev.off()


corGSVA=cor(t(geneSetEnrichSortZ))
hCor=hclust(as.dist((1-corGSVA)/2))

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Astrocytes_unannotated/GSVA_corr_v2.png",sep='')
png(pngPath, width=6,height=4,units="in",res=1200, pointsize = 10)

heatmap3(corGSVA,
         col=barColorsCor, breaks=breakBarColorsCor,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",
         labRow = colnames(astroGeneSet), labCol = NA, main="GSVA for neuron gene sets") 

dev.off()

blah=DataGSVA.astrocytes


#Permutation p-value
if (compPermute==1){
  R=1000
  meanOut3=matrix(NA,nrow=length(astroGeneSet), ncol = R)
  meanOut2=matrix(NA,nrow=length(astroGeneSet), ncol = R)
  meanOut1=matrix(NA,nrow=length(astroGeneSet), ncol = R)
  
  for (i in 1:R)
  {
    blah2=blah
    rownames(blah2)=rownames(blah)[sample(1:dim(blah)[1])]
    gsvaBoot=gsva(blah2,astroGeneSet, mx.diff=TRUE, kcdf="Gaussian")
    
    
    meanOut3[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPurple]) #mean C2 - Ad2
    meanOut2[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indRed])    #mean C2 - AD1
    meanOut1[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indBlue])   #mean C2 - C1
    
    
    print(i)
  }
  meanTrue3=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPurple]) #mean C2 - AD2
  meanTrue2=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indRed])   #mean C2 - AD1
  meanTrue1=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indBlue]) #mean C2 - C1
  
  pBoot3=matrix(NA,nrow=length(astroGeneSet), ncol = 1)
  pBoot2=matrix(NA,nrow=length(astroGeneSet), ncol = 1)
  pBoot1=matrix(NA,nrow=length(astroGeneSet), ncol = 1)
  
  for (j in 1:length(astroGeneSet))
  {
    pBoot3[j] = mean(abs(meanOut3[j,]) > abs(meanTrue3[j])) # AD2 - C2
    pBoot2[j] = mean(abs(meanOut2[j,]) > abs(meanTrue2[j])) # AD1 - C2
    pBoot1[j] = mean(abs(meanOut1[j,]) > abs(meanTrue1[j])) # C1 - C2
    
  }
  pBootFDR3=p.adjust(pBoot3,method = "fdr",)
  pBootFDR2=p.adjust(pBoot2,method = "fdr",)
  pBootFDR1=p.adjust(pBoot1,method = "fdr",)
  
  
  pBoot=data.frame(pBoot1, pBoot2, pBoot3)
  pBootFDR=data.frame(pBootFDR1, pBootFDR2, pBootFDR3)
  
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
heatmap3(DataGSVA.neuron, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA, main="astrocyte genes") 

#Compute p-values and FDR q-values for each gene set.
numGenes=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)

for (i in 1:dim(geneSetEnrich)[1])
{
  
  indGetGS=which(match(names(astroGeneSet),row.names(geneSetEnrich)[i])==1)
  numGenes[i]=length(unlist(astroGeneSet[indGetGS]))
  
}

geneSets_Stats=data.frame(rownames(geneSetEnrichSort),pBootNeu,pBootFDRNeu, numGenes);
write.xlsx(geneSets_Stats, file=paste("geneSetStats_Astrocytes_FXN_",geneSetNumber,"_norm_unannotated.xlsx",sep=''), sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)



#Plot bar graphs and gene set heatmaps

meanBluevsOR=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)
meanBluevsOR4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)
meanBluevsBlue4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)

for (indGeneSetPlot in 1:dim(geneSetEnrich)[1])
  
{
  dataBar=data.frame(as.data.frame(geneSetEnrich[indGeneSetPlot,indSort]), as.data.frame(allColorGroups))
  colnames(dataBar)=c("ES","colors")
  
  #Summarize stats for each group color
  statsOut = summarySE(dataBar, measurevar="ES", groupvars="colors") #count/mean/SD/CI 
  
  meanBluevsOR[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered"] # C2 - Ad1
  meanBluevsOR4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered4"] # C2 - AD2
  meanBluevsBlue4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="skyblue4"] # C2 - C1
  
  #Re-order the colors to match the heatmaps
  statsOut$colors = factor(statsOut$colors, levels = c("skyblue4","skyblue","orangered","orangered4"))
  
  group.colors <- c(skyblue4 = "skyblue4", skyblue = "skyblue", orangered ="orangered", orangered4 = "orangered4")
  
  pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Astrocytes_unannotated/geneSets",indGeneSetPlot,".png",sep='') 
  png(pngPath, width=4,height=4,units="in",res=600)
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
            axis.title.x = element_text(face = "plain", color = "black", 
                                        size = 24),
            axis.text.y = element_text(face = "plain", color = "black", 
                                       size = 24),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5))+
      xlab("")+
      ylab("Enrichment Score")+
      ggtitle(rownames(geneSetEnrich)[indGeneSetPlot])+
      scale_fill_manual(values=group.colors)+
      ylim(-0.5,0.5)+
      scale_x_discrete(labels=c("C1", "C2", "A1", "A2"))+
      geom_hline(yintercept=0)
    
  })
  dev.off()
  
}
#Make ES difference bar plot for neurons

### C2 v AD1
geneSetDifferenceData.C2vA1=data.frame(meanBluevsOR,rownames(geneSetEnrich),pBootFDRNeu[,2])
colnames(geneSetDifferenceData.C2vA1)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Astrocytes_C2vAD1_rev/ES_diff_FXN_Neu_C2vAD1_",geneSetNumber,".pdf",sep='')
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
    geom_text(aes(label=format((pBootFDRNeu[,2]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
    geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vA1)[1]))
  
})
dev.off()

### C2 v AD2
geneSetDifferenceData.C2vA2=data.frame(meanBluevsOR4,rownames(geneSetEnrich),pBootFDRNeu[,3])
colnames(geneSetDifferenceData.C2vA2)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Astrocytes_C2vAD2_rev/ES_diff_FXN_Neu_C2vAD2_",geneSetNumber,".pdf",sep='')
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
    #scale_fill_manual(name = "area", values=c("red","grey50"))+
    scale_fill_manual(name = "area", values=c("red", "grey50"))+
    xlab("Gene Set")+
    ylab("ES Diff. (C2-A2)")+
    geom_hline(yintercept=0)+
    geom_text(aes(label=format((pBootFDRNeu[,3]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
    geom_text(aes(label="qFDR", y=0.8, x=16.5), size=6)+
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vA2)[1]))
  
})
dev.off()

### C2 v C1
geneSetDifferenceData.C2vC1=data.frame(meanBluevsBlue4,rownames(geneSetEnrich),pBootFDRNeu[,1])
colnames(geneSetDifferenceData.C2vC1)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Astrocytes_C2vC1_rev/ES_diff_FXN_Neu_C2vC1_",geneSetNumber,".pdf",sep='')
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
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]))
  
})
dev.off()    





#____MICROGLIA___#### @@@@ UNANNOTADED
ySortZMicroglia_ann <- ySortZMicroglia
rownames(ySortZMicroglia_ann) <- microGeneNames

#seperate data based on clusters
microDataClust.1 <- ySortZMicroglia_ann[clust.microglia.cut == 1,]
microDataClust.2 <- ySortZMicroglia_ann[clust.microglia.cut == 2,]
microDataClust.3 <- ySortZMicroglia_ann[clust.microglia.cut == 3,]
microDataClust.4 <- ySortZMicroglia_ann[clust.microglia.cut == 4,]
microDataClust.5 <- ySortZMicroglia_ann[clust.microglia.cut == 5,]

#export micrglia data with cluster annotation
writeMicroData <- as.data.frame(ySortZMicroglia_ann)
writeMicroData$symbol <- rownames(writeMicroData)
writeMicroData$cluster <- "na"

writeMicroData[clust.microglia.cut == 1,]$cluster <- "cluster1"
writeMicroData[clust.microglia.cut == 2,]$cluster <- "cluster2"
writeMicroData[clust.microglia.cut == 3,]$cluster <- "cluster3"
writeMicroData[clust.microglia.cut == 4,]$cluster <- "cluster4"
writeMicroData[clust.microglia.cut == 5,]$cluster <- "cluster5"

writeMicroData <- (writeMicroData[,c("symbol","cluster")])

write.xlsx(writeMicroData, file=paste("./unannotated/MicrogliaClusterGeneSet.xlsx",sep=''), sheetName="Sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)

#save gene names of clusters
microGenesClust.1 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 1,])
microGenesClust.2 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 2,])
microGenesClust.3 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 3,])
microGenesClust.4 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 4,])
microGenesClust.5 <- rownames(ySortZMicroglia_ann[clust.microglia.cut == 5,])


# ### GO analysis to identify Biological Functions of the clusters ###
# 
#convert gene symbols to EntrezIDs
entrezClust1 <- mapIds(org.Hs.eg.db, microGenesClust.1, "ENTREZID", "SYMBOL")
entrezClust2 <- mapIds(org.Hs.eg.db, microGenesClust.2, "ENTREZID", "SYMBOL")
entrezClust3 <- mapIds(org.Hs.eg.db, microGenesClust.3, "ENTREZID", "SYMBOL")
entrezClust4 <- mapIds(org.Hs.eg.db, microGenesClust.4, "ENTREZID", "SYMBOL")
entrezClust5 <- mapIds(org.Hs.eg.db, microGenesClust.5, "ENTREZID", "SYMBOL")
# 
# #perfrom GO 

GOclust1 <- enrichGO(entrezClust1, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE) # M1
GOclust2 <- enrichGO(entrezClust2, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE) # M2
GOclust3 <- enrichGO(entrezClust3, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE) # M5
GOclust4 <- enrichGO(entrezClust4, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE) # M3
GOclust5 <- enrichGO(entrezClust5, 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.01 , pAdjustMethod = "fdr", readable= TRUE) # M4



#plot GO results 
png("./GO/GO_micro_M3.png")
barplot(GOclust1)
dev.off()

png("./GO/GO_micro_M4.png")
barplot(GOclust2)
dev.off()

png("./GO/GO_micro_M1.png")
barplot(GOclust3)
dev.off()

png("./GO/GO_micro_M2.png")
barplot(GOclust4)
dev.off()

png("./GO/GO_micro_M5.png")
barplot(GOclust5)
dev.off()

# A - clust 3 (405) - M1 
# B - clust 4 (257) - M2
# C - clust 1 (260) - M3
# D - clust 2 (277) - M4
# E - clust 5 (196) - M5

# 
# rosmap_micro_genes <- read.xlsx("../ClusterROSMAP/MicrogliaClusterGeneSet.xlsx", sheetIndex = 1)
# sinai_micro_genes <- read.xlsx("../ClusterMtSinai/MicrogliaClusterGeneSet.xlsx", sheetIndex = 1)
# 
# 
# common <- Reduce(intersect, list(microGeneNames, rosmap_micro_genes$symbol),sinai_micro_genes$symbol)
# write.xlsx(common, file=paste("./CommonMicrogliaGenes.xlsx",sep=''), sheetName="Sheet1", col.names=FALSE, row.names=FALSE, append=FALSE)
# 
# edox1 <- pairwise_termsim(GOclust1, nClusters=2)
# treeplot(BPtree1)
# 
# edox1 <- pairwise_termsim(GOclust1)
# treeplot(edox1)
# 
# edox1 <- pairwise_termsim(GOclust1)
# treeplot(edox1)
# 
# edox1 <- pairwise_termsim(GOclust1)
# treeplot(edox1)
# 
# 
# edox <- pairwise_termsim(enrichall, method = "JC")
# tree <- treeplot(edox, hclust_method = "ward.D2")
# edox <- pairwise_termsim(enrichall, method = "average")


# A - clust 3 (405) - M1 
# B - clust 4 (257) - M2
# C - clust 1 (260) - M3
# D - clust 2 (277) - M4
# E - clust 5 (196) - M5



#generate geneSet for GSVA
microGenes.total <- list(microGenesClust.1,microGenesClust.2,microGenesClust.3,microGenesClust.4,microGenesClust.5)
list.names <- c("M-3","M-4","M-1","M-2","M-5")
names(microGenes.total) <- list.names
microGeneSet =lapply(microGenes.total, function(x) x[!is.na(x)])

# DATA for GSVA in Microglia 
#DataGSVA.microglia <- yMicrogliaClean #are the normalized data (after vsn) but only of microglia (clust = 17) // NOT z-scored
DataGSVA.microglia <- ySortNormMicrogliaGenes #are the normalized data (after vsn) but only of microglia (clust = 17) // NOT z-scored
rownames(DataGSVA.microglia)  <- microGeneNames

geneSetEnrich=gsva(DataGSVA.microglia, microGeneSet, mx.diff=TRUE, kcdf="Gaussian")

geneSetEnrich_Microglia=geneSetEnrich


#----GSVA heatmap----

geneSetEnrichSort=geneSetEnrich_Microglia[,indSort] #sort samples by phenotype
geneSetEnrichSortZ=t(apply(geneSetEnrichSort,1,scale)); #z-score the data

hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Microglia/GSVA_cluster_v2.png",sep='')
png(pngPath, width=5,height=4,units="in",res=1200, pointsize = 16)

print({
  heatmap3(geneSetEnrichSortZ, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
           col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=4.5), 
           Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
           labRow = colnames(microGeneSet), labCol = NA, main="GSVA for microglia gene sets") 
})
dev.off()


corGSVA=cor(t(geneSetEnrichSortZ))
hCor=hclust(as.dist((1-corGSVA)/2))

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Microglia/GSVA_corr_v2.png",sep='')
png(pngPath, width=6,height=4,units="in",res=1200, pointsize = 10)
print({
  heatmap3(corGSVA,
           col=barColorsCor, breaks=breakBarColorsCor,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
           Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",
           labRow = colnames(microGeneSet), labCol = NA, main="GSVA for Astrocyte gene sets") 
})
dev.off()



blah=DataGSVA.microglia


if (compPermute==1){
  R=1000
  meanOut3=matrix(NA,nrow=length(microGeneSet), ncol = R)
  meanOut2=matrix(NA,nrow=length(microGeneSet), ncol = R)
  meanOut1=matrix(NA,nrow=length(microGeneSet), ncol = R)
  
  for (i in 1:R)
  {
    blah2=blah
    rownames(blah2)=rownames(blah)[sample(1:dim(blah)[1])]
    gsvaBoot=gsva(blah2,microGeneSet, mx.diff=TRUE, kcdf="Gaussian")
    
    meanOut3[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indPurple]) #mean C2 - Ad2
    meanOut2[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indRed])    #mean C2 - AD1
    meanOut1[,i]=rowMeans(gsvaBoot[,indGreen])-rowMeans(gsvaBoot[,indBlue])   #mean C2 - C1
    
    
    print(i)
  }
  meanTrue3=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indPurple]) #mean C2 - AD2
  meanTrue2=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indRed])   #mean C2 - AD1
  meanTrue1=rowMeans(geneSetEnrich[,indGreen])-rowMeans(geneSetEnrich[,indBlue]) #mean C2 - C1
  
  pBoot3=matrix(NA,nrow=length(microGeneSet), ncol = 1)
  pBoot2=matrix(NA,nrow=length(microGeneSet), ncol = 1)
  pBoot1=matrix(NA,nrow=length(microGeneSet), ncol = 1)
  
  for (j in 1:length(microGeneSet))
  {
    pBoot3[j] = mean(abs(meanOut3[j,]) > abs(meanTrue3[j])) # AD2 - C2
    pBoot2[j] = mean(abs(meanOut2[j,]) > abs(meanTrue2[j])) # AD1 - C2
    pBoot1[j] = mean(abs(meanOut1[j,]) > abs(meanTrue1[j])) # C1 - C2
    
  }
  pBootFDR3=p.adjust(pBoot3,method = "fdr",)
  pBootFDR2=p.adjust(pBoot2,method = "fdr",)
  pBootFDR1=p.adjust(pBoot1,method = "fdr",)
  
  
  pBoot=data.frame(pBoot1, pBoot2, pBoot3)
  pBootFDR=data.frame(pBootFDR1, pBootFDR2, pBootFDR3)
  
  rownames(pBoot)=rownames(geneSetEnrich)
  rownames(pBootFDR)=rownames(geneSetEnrich)
  #End permutation
  
  pBootMic=pBoot
  pBootFDRMic=pBootFDR
}

if(geneSetNumber=="1"){
  save_pBootAst1=data.frame(pBootMic,pBootFDRMic)
}else{
  save_pBootAst2=data.frame(pBootMic,pBootFDRMic)
}


#plotting astrocyte genes
heatmap3(DataGSVA.microglia, ColSideColors =combinedColorALL[indSort,], ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA, main = "Microglia Genes") 


numGenes=matrix(, nrow = dim(geneSetEnrichSort)[1], ncol = 1)

for (i in 1:dim(geneSetEnrich)[1])
{ 
  
  indGetGS=which(match(names(microGeneSet),row.names(geneSetEnrich)[i])==1)
  numGenes[i]=length(unlist(microGeneSet[indGetGS]))
  
}


geneSets_Stats=data.frame(rownames(geneSetEnrichSort), pBootMic, pBootFDRMic, numGenes);
write.xlsx(geneSets_Stats, file=paste("geneSetStats_Microglia_FXN_",geneSetNumber,"_norm.xlsx",sep=''), sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)



#Plot bar graphs and gene set heatmaps

meanBluevsOR=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)
meanBluevsOR4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)
meanBluevsBlue4=matrix(, nrow = dim(geneSetEnrich)[1], ncol = 1)

for (indGeneSetPlot in 1:dim(geneSetEnrich)[1])
  
{
  dataBar=data.frame(as.data.frame(geneSetEnrich[indGeneSetPlot,indSort]), as.data.frame(allColorGroups))
  colnames(dataBar)=c("ES","colors")
  
  #Summarize stats for each group color
  statsOut = summarySE(dataBar, measurevar="ES", groupvars="colors")
  
  
  meanBluevsOR[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered"] # AD1 - C2
  meanBluevsOR4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="orangered4"] # AD2 - C2
  meanBluevsBlue4[indGeneSetPlot]=statsOut$ES[statsOut$color=="skyblue"]-statsOut$ES[statsOut$color=="skyblue4"] # C1 - C2
  
  #Re-order the colors to match the heatmaps
  statsOut$colors = factor(statsOut$colors, levels = c("skyblue4","skyblue","orangered","orangered4"))
  
  group.colors <- c(skyblue4 = "skyblue4", skyblue = "skyblue", orangered ="orangered", orangered4 = "orangered4")
  
  pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Microglia_unannotated/geneSets_",indGeneSetPlot,".png",sep='')
  png(pngPath, width=4,height=4,units="in",res=600)
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
            axis.title.x = element_text(face = "plain", color = "black",
                                        size = 24),
            axis.text.y = element_text(face = "plain", color = "black",
                                       size = 24),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5))+
      xlab("")+
      ylab("Enrichment Score")+
      ggtitle(rownames(geneSetEnrich)[indGeneSetPlot])+
      scale_fill_manual(values=group.colors)+
      ylim(-1,1)+
      scale_x_discrete(labels=c("C1", "C2", "A1", "A2"))+
      geom_hline(yintercept=0)
    
  })
  dev.off()

}

#Make ES difference bar plot

### C2 v AD1
geneSetDifferenceData.C2vA1=data.frame(meanBluevsOR,rownames(geneSetEnrich),pBootFDRMic[,2])
colnames(geneSetDifferenceData.C2vA1)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Microglia_C2vAD1_rev/ES_diff_FXN_Ast_C2vAD1_",geneSetNumber,".pdf",sep='')
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
    geom_text(aes(label=format((pBootFDRMic[,2]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
    geom_text(aes(label="qFDR", y=0.8, x=18.5), size=6)+
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vA1)[1]))
  
})
dev.off()

### C2 v AD2
geneSetDifferenceData.C2vA2=data.frame(meanBluevsOR4,rownames(geneSetEnrich),pBootFDRMic[,3])
colnames(geneSetDifferenceData.C2vA2)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Microglia_C2vAD2_rev/ES_diff_FXN_Ast_C2vAD2_",geneSetNumber,".pdf",sep='')
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
    geom_text(aes(label=format((pBootFDRMic[,3]), scientific=FALSE, digits=4), y=.8), vjust=0.3) +
    geom_text(aes(label="qFDR", y=0.8, x=18.5), size=6)+
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vA2)[1]))
  
})
dev.off()

### C2 v C1
geneSetDifferenceData.C2vC1=data.frame(meanBluevsBlue4,rownames(geneSetEnrich),pBootFDRMic[,1])
colnames(geneSetDifferenceData.C2vC1)=c("difference","geneSetName", 'pBootFDR')

pngPath=paste("./unannotated/FXN",geneSetNumber,"_norm_Microglia_C2vC1_rev/ES_diff_FXN_Ast_C2vC2_",geneSetNumber,".pdf",sep='')
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
    geom_text(aes(label="qFDR", y=0.8, x=18.5), size=6)+
    coord_flip(ylim = c(-1,1), xlim=c(1,dim(geneSetDifferenceData.C2vC1)[1]))
  
})
dev.off()  



