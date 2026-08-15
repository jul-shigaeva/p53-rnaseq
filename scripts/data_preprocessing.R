# скачивание и активация библиотек

library("tidyverse")

# загрузка данных 

expr <- readRDS("data/raw_counts.rds")
metadata <- readRDS("data/geo_metadata.rds")

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

metadata_rnaseq %>%
  filter(!str_detect(title, "IR 10Gy$")) %>%
  select(geo_accession, title)

unique(str_extract(
  metadata_rnaseq$title,
  "(?<=h, ).+(?=, rep)"
))

metadata_rnaseq <- metadata_rnaseq %>%
  mutate(
    treatment = str_extract(title, "(?<=h, ).+(?=, rep)"),
    treatment = str_replace(treatment, "\\s+$", "")
  )

metadata_rnaseq
table(metadata_rnaseq$treatment, useNA = "ifany")

metadata_fig1 <- metadata_rnaseq %>%
  filter(treatment == "IR 10Gy" | is.na(treatment))
