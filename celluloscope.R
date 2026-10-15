library(readxl)
library(dplyr)
library(stringr)
library(writexl)
library(tidyr)
library(purrr)
library(readr)

#Leer archivo
df <- readxl::read_excel("hmmer_minor_results+1.xlsx", col_names = TRUE)
head(df)

cat("Filas al cargar archivo: ", nrow(df), "\n")

#Renombrar las columnas
colnames(df) <- c(
  "target_name", "accession", "acc", "acc_clean", "name",
  "tlen", "query_name", "acce2", "qlen", "Evalue_full", "score_full", "bias_full",
  "instance", "of", "cEvalue_dom", "iEvalue_dom", "score_dom", "bias_dom",
  "hmm_from", "hmm_to", "ali_from", "ali_to",
  "env_from", "env_to", "acc2", "description_target",
  "DESCRIPTION_UNIPROT", "Domain"
)

#Quitar columnas vacías
df <- df[, !is.na(colnames(df)) & colnames(df) != ""]

#Asegurar las variables numericas
df <- df %>%
  mutate(
    ali_from = as.numeric(ali_from),
    ali_to   = as.numeric(ali_to)
  )

#Asignar categorías grandes
df <- df %>% mutate(domain_class = case_when(
  name == "cohesin" ~ "cohesin",
  name == "cbm" ~ "cbm",
  name == "slh" ~ "slh",
  name == "dockerin" ~ "dockerin",
  name == "alt_anchor" ~ "alt_anchor",
  name %in% c("gh", "ce", "pl") ~ "enzyme",
  name == "rdgl" ~ "regulator",
  name == "x" ~ "x_domain",
  TRUE ~ "other"
))

cat("Filas después de clasificar domain_class: ", nrow(df), "\n")

#FASE 1: Remover falsos positivos

##Definir combinaciones a descartar:
remove_patterns <- list(
  c("enzyme"),
  c("cbm"),
  c("slh"),
  c("alt_anchor"),
  c("enzyme","regulator"),
  c("enzyme","slh"),
  c("enzyme","cbm"),
  c("cbm","slh"),
  c("cbm","regulator"),
  c("enzyme","cbm","slh"),
  c("slh","regulator"),
  c("x_domain"),
  c("x_domain","enzyme"),
  c("x_domain","enzyme","cbm"),
  c("x_domain","regulator"),
  c("x_domain","slh"),
  c("alt_anchor","enzyme"),
  c("regulator"),
  c("alt_anchor","regulator"),
  c("alt_anchor","slh"),
  c("dockerin","slh"),
  c("cbm","dockerin","slh"),
  c("regulator","slh","dockerin"),
  c("x_domain","slh","dockerin"),
  c("alt_anchor","dockerin","x_domain"),
  c("cbm","regulator","cohesin"),
  c("regulator","cohesin"),
  c("regulator","dockerin"),
  c("regulator","x_domain","cohesin"),
  c("regulator","x_domain","dockerin"),
  c("regulator","slh","cohesin"),
  c("enzyme","cbm","cohesin"),
  c("x_domain"),
  c("x_domain","cbm"),
  c("alt_anchor","slh","cbm","x_domain"),
  c("alt_anchor","slh","x_domain")
)

##Definir combinaciones válidas:
keep_patterns <- list(
  c("cohesin"),
  c("dockerin"),
  c("cohesin","dockerin"),
  c("enzyme","dockerin"),
  c("enzyme","cohesin","dockerin"),
  c("enzyme","cohesin"),
  c("enzyme","cbm","dockerin"),
  c("enzyme","cbm","dockerin","x_domain"),
  c("enzyme","cbm","dockerin","cohesin","x_domain"),
  c("cohesin","cbm"),
  c("cbm","dockerin"),
  c("cohesin","x_domain"),
  c("cohesin","slh"),
  c("dockerin","x_domain"),
  c("cohesin","cbm","dockerin"),
  c("alt_anchor","cohesin"),
  c("cohesin","dockerin","x_domain"),
  c("x_domain","slh","cohesin"),
  c("x_domain","cbm","cohesin"),
  c("x_domain","cbm","dockerin"),
  c("x_domain","enzyme","dockerin"),
  c("regulator"),
  c("slh","cohesin","dockerin"),
  c("dockerin","alt_anchor")
)

##Función auxiliar para comparar exactamente
matches_exact <- function(domains, patterns_list) {
  any(map_lgl(patterns_list, ~ setequal(.x, domains)))
}

##Filtrar proteínas segun reglas exactas
domains_by_protein <- df %>%
  group_by(query_name) %>%
  summarise(domains = list(unique(domain_class)))

domains_by_protein <- domains_by_protein %>%
  mutate(
    remove = map_lgl(domains, ~ matches_exact(.x, remove_patterns)),
    keep   = map_lgl(domains, ~ matches_exact(.x, keep_patterns))
  )

proteins_to_keep <- domains_by_protein %>%
  filter(!remove | keep) %>%
  pull(query_name)

df <- df %>% filter(query_name %in% proteins_to_keep)

cat("Filas después de aplicar reglas de combinación (fase 1):", nrow(df), "\n")


#FASE 2

##Funciones para evaluar solapamiento
is_partial_overlap <- function(a1, a2, b1, b2){
  (a1 < b2) & (b1 < a2) & !(a1 >= b1 & a2 <= b2) & !(b1 >= a1 & b2 <= a2)
}

is_contained_in <- function(a1, a2, b1, b2){
  (a1 >= b1) & (a2 <= b2)
}

is_container <- function(a1, a2, b1, b2){
  (a1 <= b1) & (a2 >= b2)
}

##Detectar solapamientos dentro de cada proteína y dominio
df_overlap <- df %>%
  arrange(query_name, domain_class, ali_from) %>%
  group_by(query_name, domain_class) %>%
  mutate(id = row_number()) %>%
  ungroup()

## Crear pares (A vs B)
pairs <- df_overlap %>%
  ### preparar dataset A
  rename(
    ali_from_a = ali_from,
    ali_to_a   = ali_to,
    id_a       = id,
    target_name_a = target_name
  ) %>%
  inner_join(
    ### preparar dataset B
    df_overlap %>% rename(
      ali_from_b = ali_from,
      ali_to_b   = ali_to,
      id_b       = id,
      target_name_b = target_name
    ),
    by = c("query_name", "domain_class"),
    relationship = "many-to-many"
  ) %>%
  filter(id_a != id_b)   # excluir comparación consigo mismo


pairs <- pairs %>% mutate(
  partial = is_partial_overlap(ali_from_a, ali_to_a, ali_from_b, ali_to_b),
  contained = is_contained_in(ali_from_a, ali_to_a, ali_from_b, ali_to_b),
  container = is_container(ali_from_a, ali_to_a, ali_from_b, ali_to_b)
)

##Crear etiqueta descriptiva por dominio
pairs <- pairs %>%
  mutate(
    coord_a = paste0(target_name_a, "(", ali_from_a, "-", ali_to_a, ")"),
    coord_b = paste0(target_name_b, "(", ali_from_b, "-", ali_to_b, ")")
  )

label_table <- pairs %>%
  group_by(query_name, domain_class, id_a) %>%
  summarise(
    overlap_info = paste(c(
      if (any(partial)) paste0("partial_overlap:", paste(coord_b[partial], collapse=",")),
      if (any(contained)) paste0("contained_in:", paste(coord_b[contained], collapse=",")),
      if (any(container)) paste0("contains_other:", paste(coord_b[container], collapse=","))
    ), collapse=";"),
    .groups = "drop"
  )

df_out <- df_overlap %>%
  left_join(label_table, by = c("query_name","domain_class","id"="id_a")) %>%
  rename(overlap = overlap_info)

df_out$overlap[is.na(df_out$overlap)] <- "none"


#Exportar archivo final
write_xlsx(df_out, "Filtered_Domains_With_Overlap+Evalue_RF1.xlsx")

cat("\nArchivo generado: Filtered_Domains_With_Overlap+Evalue_RF1.xlsx\n")

