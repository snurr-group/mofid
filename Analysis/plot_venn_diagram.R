# Analyzing the overlap between GA hMOF, ToBaCCo, and CoRE databases
# Imports data from find_duplicates.sh and runs all of the necessary joins to aggregate
# duplicates results into a Venn diagram (of unique MOFs, excluding topology error codes)

library(readr)
library(dplyr)
library(tidyr)

library(VennDiagram)
library(cowplot)


DUPLICATES_PATH <- "Analysis/Figures/Duplicates/"
read_dup_data <- function(db_description) {
  path <- file.path(DUPLICATES_PATH, paste("overlap", db_description, "summary.tsv", sep="_"))
  read_tsv(path, col_types="ciicci")
}

# Overlap data (between databases)
o_tob_ga <- read_dup_data("tob_and_ga")
o_core_ga <- read_dup_data("with_ga")
o_core_tob <- read_dup_data("with_tob")

# Duplicates data (within database), so we can base our results off of "unique MOFs"
get_num_lines <- function(db_name) {
  read_tsv(
    file.path(DUPLICATES_PATH, paste("duplicates", db_name, "all_families.tsv", sep="_")),
    col_types = "cic"
    ) %>% 
    nrow()
}
num_unique_mofs <- c(
  core = get_num_lines("core"),
  ga = get_num_lines("ga"),
  tob = get_num_lines("tob")
)

# Find the MOFs in common between all three databases
o_combined <- 
  select(o_core_ga, identifier, ga_qty = right_qty) %>% 
  inner_join(select(o_core_tob, identifier, tob_qty = right_qty, tob_name = right_filenames), by="identifier")

# Construct Venn diagrams
# See resources at https://community.rstudio.com/t/how-to-create-a-venn-diagram/19228/5
# which links to a good vignette at https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
grid.newpage()
draw.triple.venn(
  area1 = num_unique_mofs["core"],
  area2 = num_unique_mofs["ga"],
  area3 = num_unique_mofs["tob"],
  n12 = nrow(o_core_ga),
  n23 = nrow(o_tob_ga),
  n13 = nrow(o_core_tob),
  n123 = nrow(o_combined),
  #category = c("CoRE MOF 2019-ASR", "GA hMOFs", "ToBaCCo"),
  category = c("CoRE MOFs", "GA hMOFs", "ToBaCCo"),
  lty = "blank",
  fill = c("skyblue", "pink1", "mediumorchid"),
  cat.dist = c(-0.05, -0.05, -0.025),
  #cat.default.pos = "text", cat.dist = c(0, 0, 0),
  cex = 1.5, cat.cex = 1.5,
  # Note: fonts in Windows are touchy.  But even if sans-serif, etc., gives warnings,
  # at least it overrides the default "serif" font.
  # Here is a reference of potentially how to do this (slowly) https://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html
  cat.fontfamily = "sans-serif", fontfamily = "sans-serif"
) %>% 
  cowplot::save_plot("Analysis/Figures/overlap_venn.png", ., base_aspect_ratio=1.5)

