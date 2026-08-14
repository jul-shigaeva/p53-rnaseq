# установка BiocManager 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# скачивание и активация библиотек

BiocManager::install("GEOquery", force = TRUE)
library("GEOquery")
library("tidyverse")
library("data.table")

# получение raw counts

url <- c("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE100099&format=file&file=GSE100099_raw_counts_GRCh38.p13_NCBI.tsv.gz")

download.file(
  url,
  destfile = "data/GSE100099_raw_counts_GRch38.p13_NCBI.tsv.gz",
  mode = "wb"
)

expr <- fread("data/GSE100099_raw_counts_GRch38.p13_NCBI.tsv.gz") # матрица каунтов 
dim(expr)
colnames(expr) # в столбцах гены, в колонках образцы (GSM)

# сохранение 

saveRDS(expr, "data/raw_counts.rds")

# получение metadata table 

gse <- getGEO("GSE100099", GSEMatrix = TRUE)
length(gse)
gse <- gse[[1]] # достаем ExpressionSet
dim(exprs(gse))
metadata <- pData(gse) # достаем таблицу с информацией об образцах 
colnames(metadata)
metadata_expr <- metadata[, c("geo_accession", "title")] # информация об образцах оказалась в колонке Title
# мы это выяснили посмотрев глазами в разные колонки 

# сохранение 
saveRDS(metadata_expr, "data/geo_metadata.rds")