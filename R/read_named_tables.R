# Reads tidy data parsed by Python/names_to_tables.py

library(dplyr)
library(readr)
library(stringr)

read_mofid <- function(folder) {
  # Ready the data in the parsed folder as a list of relevant data frames
  cat <- 
    file.path(folder, "cat.tsv") %>% 
    read_tsv(col_names=c("name", "cat"), col_types="ci", na="None")
  smiles <- 
    file.path(folder, "smiles.tsv") %>% 
    read_tsv(col_names=c("name", "smiles"), col_types="cc", na="ERROR")
  smiles_part <- 
    file.path(folder, "smiles_part.tsv") %>% 
    read_tsv(col_names=c("name", "smiles_part"), col_types="cc", na="ERROR")
  topology <- 
    file.path(folder, "topology.tsv") %>% 
    read_tsv(col_names=c("name", "topology"), col_types="cc")
  list(cat=cat, smiles=smiles, smiles_part=smiles_part, topology=topology)
}

bracket_pattern <- "\\[[A-Za-z]+\\]"

strip_core_suffix <- function(x) {
  # Greedily matches everything before the final underscore and deletes everything else
  str_replace(x, "(.*)_.*", "\\1")
}

normalize_filenames <- function(mofid_tables) {
  # Normalizes the names from the results of read_mofid to be the refcodes or DOI without cleanup info
  lapply(mofid_tables, mutate, name=strip_core_suffix(name))
}


