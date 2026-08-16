# скачивание и активация библиотек

library(tidyverse)
library(DESeq2)

# загрузка данных 

expr <- readRDS("data/expr_fig1.rds")
metadata <- readRDS("data/metadata_fig1.rds")

# приведение матрицы экспрессии к нормальному виду  

summary(expr[, -1]) # воу...

counts <- expr[, -1]
rownames(counts) <- expr$GeneID

# фильтрация низкоэкспрессируемых генов

summary(rowSums(counts > 0))

keep <- rowSums(counts >= 10) >= 2 # зафиксировать в README, неоднозначно какой критерий брать 
counts_filtered <- counts[keep, ]
dim(counts_filtered)



metadata <- metadata %>%
  mutate(time = as.numeric(time))

# нормализация 

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)

counts_norm <- counts(dds, normalized = TRUE)

# моделирование по времени 

metadata %>%
  dplyr::count(time, p53, treatment)

metadata_ordered <- metadata %>%
  arrange(time, p53, rep)
sample_order <- metadata_ordered$geo_accession
counts_ordered <- counts_norm[, sample_order]
all(colnames(counts_ordered) == metadata_ordered$geo_accession)

counts_log <- log2(counts_norm + 1)
summary(as.vector(counts_log))

counts_z <- t(scale(t(counts_log)))
dim(counts_z)


set.seed(123)

genes_test <- sample(rownames(counts_z), 100)

heatmap_test <- counts_z[genes_test, ]

library(pheatmap)

pheatmap(
  heatmap_test,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE
)
