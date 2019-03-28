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

# Extract the MC and ON from ToBaCCo MOFs
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
df_get_mc <- Vectorize(function(x, n) {return(x[n])})
tobacco_errs <- tobacco_errs %>% 
  mutate(nodes.on = df_extract_on(code.nodes)) %>% 
  mutate(nodes.on = ifelse(is.na(nodes.on), "None", nodes.on)) %>% 
  mutate(nodes.mc = df_extract_on(code.nodes, "_mc_")) %>% 
  mutate(mc1 = ifelse(nodes.mc=="Multiple", unlist(df_get_mc(code.nodes, 1)), nodes.mc)) %>% 
  mutate(mc2 = ifelse(nodes.mc=="Multiple", unlist(df_get_mc(code.nodes, 2)), NA))


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
  # Translate error codes from check_mof_composition.py into their implications
  
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


tobacco_errs <- tobacco_errs %>% 
  mutate(err_type=df_summarize_errors(errors)) %>% 
  #mutate(err_color=ifelse(err_type == "Definition mismatch", "Mismatch", "Error")) %>% 
  # More aggressively detect Systre crashes, since the Python code will hide them if multiple sources of error
  mutate(err_type=ifelse(str_detect(from_mofid, "ERROR"), "Crash", err_type))
ga_errs <- ga_errs %>% 
  mutate(err_type=df_summarize_errors(errors)) %>% 
  mutate(err_type=ifelse(str_detect(from_mofid, "ERROR"), "Crash", err_type))


### UNDERSTANDING ERROR FLOWS ###

# TODO: consider defining classes of nodes/linkers if we use that information

RUN_ERROR_FLOWS <- TRUE
if (RUN_ERROR_FLOWS) {
  # There are two ways we can think about diagnosing the errors in the code
  # 1. Going from errors to their cause (this block)
  # 2. Going from patterns in the heatmap (e.g. certain linkers) back to types of error
  # Both approaches could be informative, but let's start with #1 here.
  
  tobacco_errs$err_type %>% unique
  colnames(tobacco_errs)
  
  # Define common classes of ToBaCCo linkers
  tob_L_categories <- tibble(
    code.linker = 1:47,  # adding the "L_" prefix below
    L_4n_ring = logical(47),  # defaults to FALSE
    L_triple_bond = logical(47),
    L_double_bond = logical(47),
    L_N_N_double = logical(47)
  )
  tob_L_categories[c(13, 26, 30, 35, 42, 44, 45), "L_4n_ring"] <- TRUE
  tob_L_categories[c(3, 16, 19, 20, 27, 31:35, 37, 38, 41, 42, 46, 47), "L_triple_bond"] <- TRUE
  tob_L_categories[c(2, 14, 15, 29, 30, 39), "L_double_bond"] <- TRUE
  tob_L_categories[c(1, 40), "L_N_N_double"] <- TRUE
  tob_L_categories <- tob_L_categories %>% mutate(code.linker = paste0("L_", code.linker))
  
  understand_tobacco <- tobacco_errs %>%
    select(nodes.mc, nodes.on, code.linker, code.topology, cat, err_type, mc1, mc2) %>% 
    left_join(tob_L_categories, by="code.linker") %>% 
    mutate_all(funs(factor))  # set as factor so we can easily run summaries
  
  # Trivial example summary for tpt vs. stp:
  understand_tobacco %>% filter(err_type=="Definition mismatch") %>% summary()
  # Now let's dig into this data more closely
  # For the figure in the paper, recall that my goal is to clearly show why these errors occur.
  # We don't necessarily need all the detail, so an "Other" column is perfectly acceptable.
  understand_tobacco %>% ggplot(aes(err_type)) + geom_bar() + coord_flip()
  # Honestly, after the analysis below, I think a bar graph might be the clearest way to make our point,
  # followed by "see text/SI for discussion of classes of error".  We just need to reorganize the bars,
  # add percentages, and possibly color the bars by success vs. node incompatibility and other common
  # sources of error, etc.
  #
  # Then again, Sankey (or a heatmap) could also clearly show the relationship between a few causes
  # and their many effects. Ask Randy and Andrew for their opinion.
  
  # By category:
  # Most of the topological errors are sym_4_mc_1 (tetrahedral zinc).  Same with crashes
  understand_tobacco %>% filter(err_type=="topology") %>% summary()
  understand_tobacco %>% filter(err_type=="Crash") %>% summary()
  # nodes errors are mostly sym_8_mc_7.  When including results from mc1/2, it's 273/282
  understand_tobacco %>% filter(err_type=="node_single_bonds") %>% summary()
  # "Many errors" is largely L_29 (strangely a benzene with double bonds coming off it)
  # Ah, and again sym_8_mc_7 problems propagating back to other sources of error 
  understand_tobacco %>% filter(err_type=="Many errors") %>% summary()
  tobacco_errs %>% filter(err_type=="Many errors" & code.linker=="L_29")
  # TODO: Based on "two errors" and the common problems with sym_8_mc_7, let's make a new class of error
  # aggregating node and linker BO problems.
  # But again, half of the errors are sym_8_mc_7.
  understand_tobacco %>% filter(err_type=="Two errors") %>% summary
  tobacco_errs %>% filter(err_type=="Two errors")
  # Looking into the linkers, nothing particularly stands out from SB's
  understand_tobacco %>% filter(err_type=="linker_single_bonds") %>% summary
  # But 564/841 of the bond order problems are these 4-N rings
  understand_tobacco %>% filter(err_type=="linker_bond_orders") %>% summary
  # Interestingly, all of the "formula" errors are L_8,9,10, which have big rings next to a phenyl.
  understand_tobacco %>% filter(err_type=="formula") %>% summary
  
  # Now repeating for the GA MOFs
  
  ga_errs$err_type %>% unique
  colnames(ga_errs)
  understand_ga <- ga_errs %>%
    select(code.nodes, code.linker1, code.linker2, cat, err_type) %>% 
    #left_join(tob_L_categories, by="code.linker") %>% 
    mutate_all(funs(factor))  # set as factor so we can easily run summaries
  
  understand_ga %>% ggplot(aes(err_type)) + geom_bar() + coord_flip()
  
  # From the bar graph, crashes are the most important to address.
  # Interestingly, 101/129 are short L_0,1,2 secondary linkers
  understand_ga %>% filter(err_type=="Crash") %>% summary()
  # All of the definition mismatches are from Zr MOFs
  understand_ga %>% filter(err_type=="Definition mismatch") %>% summary()
  # Many of the "two error" cases are the definition mismatch plus something else
  understand_ga %>% filter(err_type=="Definition mismatch") %>% summary()
  ga_errs %>% filter(err_type=="Two errors")
  # Most of the single bonds errors are the gigantic L_4 linker
  understand_ga %>% filter(err_type=="linker_single_bonds") %>% summary()
  # 7/15 of the linker_bond_orders errors are the 4N linker L_27
  understand_ga %>% filter(err_type=="linker_bond_orders") %>% summary()
  understand_ga %>% filter(err_type=="linker_bond_orders") %>%
    filter(code.linker1==37 | code.linker2==37) %>% nrow
  # Overall, only 260 non-successes in this DB, including the definition mismatches
  understand_ga %>% filter(err_type!="Success") %>% nrow
  
}  # endif(RUN_ERROR_FLOWS)


### PLOTS FOR MAIN TEXT ###

# TODO: implement figure for the paper, whether it's Sankey, a heatmap, or a simple filled bar chart

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

