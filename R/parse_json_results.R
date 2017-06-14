# Parse JSON results files from Python/check_mof_linkers.py

library(jsonlite)
library(dplyr)
library(stringr)

import_err_df <- function(filename) {
  df <- fromJSON(filename)$mofs
  code_info <- df$`_codes`
  if (!is.null(code_info)) {
    df$`_codes` <- NULL
    colnames(code_info) <- paste0("code.", colnames(code_info))
    df <- bind_cols(df, code_info)
  }
  df %>% select(name, errors, everything())
}

inlist <- function(s, x) {
  # Applies the %in% command to a vector of list of element(s)
  lapply(
    x,
    function(y) {
      s %in% unlist(y)
      }
    ) %>% 
  unlist  # lapply returns a list, but we want a vector of logicals
}

length_list <- function(x) {
  # For each element in a vector of list of elements, gets the length of the contained vector
  lapply(
    x,
    function(y) {length(unlist(y))}
  ) %>% 
  unlist
}

get_before_space <- function(x) {
  # What is the content before a space in a column?
  # E.g. if you have a column of full SMILES, what is the chemical part of that column?
  # Example: three %>% filter(inlist("err_smiles", errors)) %>% mutate(cif_smiles = get_before_space(from_cif)) %>% View
  lapply(
    x,
    function(y) {
      str_split(y, " ")[[1]][1]
      }
    ) %>%
    unlist
}


# Examples:
# filter(import_err_df("Notebooks/20170605-tobacco/very_first_big_tobacco.json"), inlist("err_smiles", errors)) %>% View
# "Notebooks/20170605-tobacco/tobacco_with_codes.json" %>% import_err_df %>% filter(inlist("sym_8_mc_7", code.nodes)) %>% View
# Could also consider using `length` to detect `character(0)`, etc.

