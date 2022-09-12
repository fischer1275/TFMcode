setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#install.packages("xlsx")
library("xlsx")


#file = system.file("tests", "Zhang2016GeneExpressionCluster_NO_NORM3.xlsx", package = "xlsx")


res=read.xlsx("EnsembleIDs.xlsx", 1, header=TRUE)

#head(res[, 1])

source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Biobase)
library(EnsDb.Hsapiens.v75)

library(EnsDb.Hsapiens.v86)


#names(res) = res$NA.;
#rownames(res) = res$ensembl_id;

#res2=as.matrix(res$ensembl_id);

#res3=as.character(res2)


res4<- transform(res, ensembl_id = as.character(ensembl_id))

entrez2 <- AnnotationDbi::select(org.Hs.eg.db, keys=res4$ensembl_id, columns='SYMBOL', keytype='ENSEMBL', multiVals='LIST')
entrez2$dups=duplicated(entrez2)

BCLs <- genes(res4$ensembl_id, columns = c("gene_name", "entrezid", "gene_biotype"),return.type = "DataFrame")

              


entrez2 <- select(EnsDb.Hsapiens.v86, keys=res4$ensembl_id, columns='SYMBOL', keytype='GENEID', multiVals='list')
entrez2$dups=duplicated(entrez2$GENEID)

write.xlsx(entrez2, file="entrezout4.xlsx", sheetName="EntrezOut", col.names=TRUE, row.names=TRUE, append=FALSE)
