# установка BiocManager 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# скачивание и активация библиотек

BiocManager::install("GEOquery", force = TRUE)
library("GEOquery")
library("tidyverse")

# загрузка данных 

expr <- readRDS("data/raw_counts.rds")
metadata <- readRDS("data/geo_metadata.rds")

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