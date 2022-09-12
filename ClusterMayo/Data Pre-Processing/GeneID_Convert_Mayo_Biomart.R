#install.packages("xlsx")

#---- preliminaries ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("xlsx")

source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
library(AnnotationDbi)
library(Biobase)
library(biomaRt)

#---- Conversion ----

res=read.xlsx("EnsembleIDs.xlsx", 1, header=TRUE)

mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
res4= transform(res, ensembl_id = as.character(ensembl_id))
G_list = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=res4,mart= mart, uniqueRows=FALSE)
write.xlsx(G_list, file="symbolOut.xlsx", sheetName="SymbolOut", col.names=TRUE, row.names=TRUE, append=FALSE)
