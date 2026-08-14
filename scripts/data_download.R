# установка BiocManager 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# скачивание библиотеки GEOquery

BiocManager::install("GEOquery", force = TRUE)
library("GEOquery")

# получение raw counts

library(data.table)
expr <- fread("C:/Users/Юлия/Downloads/GSE100099_raw_counts_GRCh38.p13_NCBI.tsv.gz")
dim(expr)
colnames(expr)

# получение metadata table 


gse <- getGEO("GSE100099", GSEMatrix = TRUE)
length(gse)
gse <- gse[[1]] # достали ExpressionSet
dim(exprs(gse))
metadata <- pData(gse) # достали таблицу с информацией об образцах 
colnames(metadata)
