# Sanity check to see how many isotherms NIST has for the CoRE MOFs

library(jsonlite)
library(dplyr)
library(readr)
library(stringr)

isotherm_list <- fromJSON("Data/NISTcheck/isotherms.json")
core_dois <- read_tsv("Data/NISTcheck/structure-doi.tsv",
                      na = c("", "NA", "-")
                      )

useful_isotherms <- inner_join(core_dois, isotherm_list)

# Alternatively, let's consider json from the bibliography
biblio_list <- fromJSON("Data/NISTcheck/biblio.json")
useful_biblio <- inner_join(core_dois, biblio_list) %>%
  select(REFCODE, DOI, category, adsorbentMaterial, adsorbateGas, temperature, pressure)
useful_expt <- useful_biblio[str_detect(useful_biblio$category, "exp"),]

# Overall, it looks like there's 215 CoRE MOFs with expt data in SRD-205
# Possibly more if you consider duplicates by name and/or structure
