# Visualize sources of error from round-tripping

# Comparing MOFid's of GA/ToBaCCo MOFs against their recipes to check for
# systematic sources of error in MOFid.
# See README for usage and prerequisites.

# In a previous talk about this project (and a notebook from 2017-08-31), I reported
# the percent success rate for the validator according to the following classes:
# - ToBaCCo MOFs: successes.  We get an additional 2.6% up to 78.2% if we accept the stp/tpt bug.
# - GA MOFs: success, replaced_pillarX, or unk_pillarX.  That means the MOF is either an exact match or
#   ambiguous results from primary vs. secondary linker for carboxylate vs. pyridine pillar.
#
# Thinking about these results for a new analysis, generate two types of figures:
# 1. For the main text, draw a Sankey/alluvial diagram showing flow of BB's vs. error types
#    See also https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
# 2. For the SI, consider using heatmaps, and possibly a sum over different error types
#    This might not be necessary -- TBD.
#
# If I need treemaps later, consider library(treemapify) and this layout:
# https://support.office.com/en-us/article/create-a-treemap-chart-in-office-dfe86d28-a610-4ef5-9b30-362d5c624b68


source('Analysis/R-Utilities/parse_json_results.R')

library(ggplot2)
library(stringr)
library(cowplot)
library(viridis)
library(ggalluvial)
library(purrr)


### DATA IMPORT ###

# Import data
tobacco_errs <- import_err_df("Summary/tob_validation.json") %>% 
  # rename the linker column removing the L_
  mutate(linker = ifelse(str_detect(code.linker, "^L_"), str_sub(code.linker, 3), code.linker))
ga_errs <- import_err_df("Summary/ga_validation.json")

# Extract the ON from ToBaCCo MOFs
extract_on <- function(list_of_nodes, pattern="_on_") {
  matches <- str_detect(list_of_nodes, pattern)
  if (sum(matches) == 0) {
    return(NA_character_)
  } else if (sum(matches) > 1) {
    return("Multiple")
  } else {
    return(list_of_nodes[matches])
  }
}
# See https://community.rstudio.com/t/dplyr-alternatives-to-rowwise/8071 on how to mutate by row
df_extract_on <- Vectorize(extract_on)
tobacco_errs <- tobacco_errs %>% 
  mutate(nodes.on = df_extract_on(code.nodes)) %>% 
  mutate(nodes.on = ifelse(is.na(nodes.on), "None", nodes.on))


### INTERPRET ERRORS ###

# Extract more detailed error information
df_nested_c_length <- Vectorize(function(column) {length(column)})
# TODO: consider lumping as part of a bigger error interpretation function
tobacco_errs <- mutate(tobacco_errs, err_qty = df_nested_c_length(errors))

# First, what are the unique error classes?
err_classes_tob <- flatten_chr(tobacco_errs$errors) %>% unique
err_classes_ga <- flatten_chr(ga_errs$errors) %>% unique
# c(err_classes_ga, err_classes_tob) %>% unique

mismatch_codes <- c("V_incomplete_linker",
  "Zr_mof_as_hex", "Zr_mof_not_fcu", "Zr_mof_pcu_as_fcu",
  "stp_from_tpt"
)
pillar_codes <- c("unk_pillar1", "unk_pillar2", "replaced_pillar1", "replaced_pillar2")

summarize_errors <- function(error_char_vec) {
  # Translate error codes from check_mof_linkers.py into their implications
  
  # First handle the easiest cases
  if (length(error_char_vec) == 0) {
    return("Success")
  }
  
  # Filter out known issues in structure recipes
  cleaned_errs <- error_char_vec[!(error_char_vec %in% pillar_codes)]
  if (length(cleaned_errs) == 0) {
    # Ambiguous pillars should not be considered a source of error
    # return("Ambiguous pillars")
    return("Success")
  }
  
  cleaned_errs <- cleaned_errs[!(cleaned_errs %in% mismatch_codes)]
  if (length(cleaned_errs) == 0) {
    return("Definition mismatch")
  }
  
  # Strip off the "_extra" suffix.  It's possibly too much granularity, and we don't need to differentiate
  # between inconsistent bond assignment and incorrect bond assignment (or treat them separately)
  cleaned_errs <- str_replace(cleaned_errs, "_extra", "") %>% unique
  
  # TODO STUFF
  if (length(cleaned_errs) == 1) {
    if (cleaned_errs %in% c("err_systre_error", "err_cpp_error")) {
      return("Crash")
    } else {
      return(substr(cleaned_errs[1], 5, nchar(cleaned_errs[1])))  # strip the "err_"
    }
  } else if (length(cleaned_errs) > 2) {
    return("Many errors")
  }
  
  # Cases when length(error_char_vec) == 2:
  return("Two errors")
}
df_summarize_errors <- Vectorize(summarize_errors)

# TODO: make a new .df to translate summarize_errors to an overall error class
# "structure_mismatch" label?
# Or just have an else statement for Success, etc.

tobacco_errs <- tobacco_errs %>% 
  mutate(err_type=df_summarize_errors(errors)) %>% 
  mutate(err_color=ifelse(err_type %in% mismatch_codes, "Mismatch", "Error"))
ga_errs <- ga_errs %>% 
  mutate(err_type=df_summarize_errors(errors)) %>% 
  mutate(err_color=ifelse(err_type %in% mismatch_codes, "Mismatch", "Error"))


### UNDERSTANDING ERROR FLOWS ###

RUN_ERROR_FLOWS <- TRUE
if (RUN_ERROR_FLOWS) {
  # TODO: run `summary` as described below to understand error flows, if I go that route
  # TO GET STARTED, take each error class and run a summary command on the relevant classes of error
}  # endif(RUN_ERROR_FLOWS)


### PLOTS FOR MAIN TEXT ###

# TODO: implement figure for the paper

# Initial overall structure for the diagram
sankey_tobacco <- 
  tobacco_errs %>%
  filter(err_type != "Success") %>%
  group_by(err_type, nodes.on, err_qty) %>%
  summarize(err_n = n()) %>%
  #ggplot(aes(y=1, axis1=linker, axis2=nodes.on, axis3=err_qty)) +
  ggplot(aes(y=err_n, axis1=err_type, axis2=nodes.on, axis3=err_qty)) +
  geom_alluvium(aes(fill=err_type)) +
  geom_stratum(width=1/6, color="grey") +
  guides(fill=FALSE) +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Error class", "on", "Num errs"), expand = c(.05, .05))
cowplot::save_plot("Analysis/Figures/sankey_tobacco.png", sankey_tobacco, base_aspect_ratio=2.0)

sankey_ga <-
  ga_errs %>%
  filter(err_type != "Success") %>%  # Report % success as a number
  group_by(err_type, code.nodes, code.linker1) %>%
  summarize(err_n = n()) %>%
  #ggplot(aes(y=1, axis1=linker, axis2=nodes.on, axis3=err_qty)) +
  ggplot(aes(y=err_n, axis1=err_type, axis2=code.nodes, axis3=code.linker1)) +
  geom_alluvium(aes(fill=err_type)) +
  geom_stratum(width=1/6, color="grey") +
  guides(fill=FALSE) +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Error class", "Metal node", "L1"), expand = c(.05, .05))
cowplot::save_plot("Analysis/Figures/sankey_ga.png", sankey_ga, base_aspect_ratio=2.0)

# TODO: consider the networkD3::sankeyNetwork layout, instead
# https://www.r-graph-gallery.com/323-sankey-diagram-with-the-networkd3-library/
# Would need to revamp the figure by running a series of steps to convert the data frame, if we went that route
# Maybe also river plot as an alternative: https://stackoverflow.com/questions/9968433/sankey-diagrams-in-r
# That package seems a bit cleaner if you don't need interactivity
# At any rate, I'll need to semi-manually calculate flows/my story instead of relying on ggplot.
# The format for makeRiver is a dataframe with N1, N2, Value.

# Idea: the purpose of the Sankey diagram will be error attribution, to complement the heatmaps
# TODO: rephrase the "success" of the heatmaps for the GA MOFs - and notate properly wihtin the text
# Also, maybe then it could be a 2-part figure for both GA/ToBaCCo.  A heatmap with overall sucess, then Sankey to drill down.


### PLOTS FOR SI ###

# TODO: consider breaking up the diagram into more specific classes of error?
# Honestly, these plots might become irrelevant if the Sankey diagram has sufficient information.

tobacco_heatmap <- 
  tobacco_errs %>% 
  group_by(linker, nodes.on) %>% 
  summarize(num_successes = sum(match), num_total = n()) %>%
  mutate(success_rate = num_successes / num_total) %>% 
  mutate(Linker = factor(linker, levels=c("_", 1:50))) %>% 
  ggplot(aes(Linker, nodes.on)) +
  geom_tile(aes(fill=success_rate)) +
  # scale_fill_viridis() +
  guides(fill=FALSE) +
  scale_x_discrete(breaks=c(1, seq(5, 47, 5)), drop=FALSE)
cowplot::save_plot("Analysis/Figures/si_recipe_tobacco.png", tobacco_heatmap)

# And repeat for catenated, non-functionalized ("parent") GA hMOFs
ga_heatmap <- 
  ga_errs %>% 
  group_by(code.linker1, code.linker2, code.nodes) %>% 
  summarize(num_successes = sum(match), num_total = n()) %>%
  mutate(success_rate = num_successes / num_total) %>% 
  mutate(Linker1 = factor(code.linker1, levels=0:39)) %>% 
  mutate(Linker2 = factor(code.linker2, levels=0:39)) %>% 
  ggplot(aes(Linker1, Linker2)) +
  geom_tile(aes(fill=success_rate)) +
  guides(fill=FALSE) +
  facet_grid(code.nodes~.) +
  scale_x_discrete(drop=FALSE) +
  scale_y_discrete(drop=FALSE)
cowplot::save_plot("Analysis/Figures/si_recipe_ga.png", ga_heatmap)

