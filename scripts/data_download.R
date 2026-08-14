# установка BiocManager 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# скачивание и активация библиотек

BiocManager::install("GEOquery", force = TRUE)
library("GEOquery")
library("tidyverse")

# получение raw counts

library(data.table)
expr <- fread("C:/Users/Юлия/Downloads/GSE100099_raw_counts_GRCh38.p13_NCBI.tsv.gz") # матрица каунтов 
dim(expr)
colnames(expr) # в столбцах гены, в колонках образцы (GSM)

# получение metadata table 

gse <- getGEO("GSE100099", GSEMatrix = TRUE)
length(gse)
gse <- gse[[1]] # достаем ExpressionSet
dim(exprs(gse))
metadata <- pData(gse) # достаем таблицу с информацией об образцах 
colnames(metadata)
metadata_expr <- metadata[, c("geo_accession", "title")] # информация об образцах оказалась в колонке Title
# мы это выяснили посмотрев глазами в разные колонки 

# проверка соответствия GSM 

gsm_expr <- colnames(expr)
gsm_metadata <- metadata$geo_accession
sum(gsm_expr %in% gsm_metadata) # на один меньше чем ожидалось 
# проверка, что это 
setdiff(colnames(expr), metadata$geo_accession)
# удаление лишней строки 
gsm_expr <- colnames(expr)[-1]
# проверка что после удаления все совпадает 
all(gsm_expr %in% metadata$geo_accession) # TRUE 

# фильтрация GSM на образцы RNA-seq 

metadata_expr %>%
  filter(str_starts(title, "RNA-Seq"))
nrow(metadata_rnaseq)

# получение информации об образцах в удобном формате 

metadata_rnaseq <- metadata_rnaseq %>%
  separate_wider_delim(
    title,
    delim = ", ",
    names = c("seq_type", "time", "treatment", "replicate"),
    too_few = "align_start" # если частей меньше, выравнивать по началу (NA справа)
  ) 

table(metadata_rnaseq$treatment, useNA = "ifany")

