# скачивание и активация библиотек

library("tidyverse")

# загрузка данных 

expr <- readRDS("data/raw_counts.rds")
metadata <- readRDS("data/geo_metadata.rds")

# РАБОТА С МЕТАДАННЫМИ 

# проверка соответствия GSM 

gsm_expr <- colnames(expr)[-1]
gsm_metadata <- metadata$geo_accession
sum(gsm_expr %in% gsm_metadata) # на 15 меньше чем ожидалось 

missing_gsm <- setdiff(gsm_metadata, gsm_expr)
metadata[metadata$geo_accession %in% missing_gsm, c("geo_accession", "title")] # это данные ChIP-seq

# фильтрация GSM на образцы RNA-seq 

metadata_rnaseq <- metadata %>%
  filter(str_starts(title, "RNA-Seq"))
nrow(metadata_rnaseq)

# получение информации об образцах в удобном формате 

metadata_rnaseq <- metadata_rnaseq %>%
  mutate(
    time = str_extract(title, "(?<=t=)[0-9.]+"), 
    rep = str_extract(title, "(?<=rep).+"), 
    p53 = str_extract(rep, "(?<=_p53).+")
  )

metadata_rnaseq <- metadata_rnaseq %>% 
  mutate(
    rep = str_extract(rep, "\\d+")
  )

metadata_rnaseq$p53[is.na(metadata_rnaseq$p53)] <- "WT"

# оставим образцы, подвергавшиеся воздействию радиации 

table(str_detect(metadata_rnaseq$title, "IR 10Gy"))

# какие варианты treatment есть 

unique(str_extract(
  metadata_rnaseq$title,
  "(?<=h, ).+(?=, rep)"
))

# создание нового столбца Treatment и форматирование к единообразию 

metadata_rnaseq <- metadata_rnaseq %>%
  mutate(
    treatment = str_extract(title, "(?<=h, ).+(?=, rep)"),
    treatment = str_replace(treatment, "\\s+$", "")
  )

table(metadata_rnaseq$treatment, useNA = "ifany")

# из статьи и дизайна эксперимента видно, что IR 10Gy + Nutlin не использовался для анализа в Figure 1
# удаление этих образцов

metadata_fig1 <- metadata_rnaseq %>%
  filter(treatment == "IR 10Gy" | is.na(treatment))

# последняя чистка 

metadata_fig1 <- metadata_fig1 %>%
  select(-title) 
rownames(metadata_fig1) <- NULL

# проверка, что для всех полученных образцов есть данные экспрессии 

all(metadata_fig1$geo_accession %in% colnames(expr)) # TRUE


# ФИЛЬТРАЦИЯ ДАННЫХ ЭКСПРЕССИИ

expr_fig1 <- expr[, c("GeneID", metadata_fig1$geo_accession)]
table(colnames(expr_fig1)[-1] == metadata_fig1$geo_accession) # проверка соответствия порядка 


# сохранение 

saveRDS(metadata_fig1, "data/metadata_fig1.rds")
saveRDS(expr_fig1, "data/expr_fig1.rds")
