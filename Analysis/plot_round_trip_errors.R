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
# Thinking about these results, the figure we need to generate is a simple bar graph showing
# the rates of success and different types of failures, along with common causes.

source('Analysis/R-Utilities/parse_json_results.R')

library(ggplot2)
library(stringr)
library(cowplot)
library(viridis)
library(ggalluvial)
library(purrr)
library(forcats)
library(RColorBrewer)
library(assertthat)


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
  mutate(mc2 = ifelse(nodes.mc=="Multiple", unlist(df_get_mc(code.nodes, 2)), "None"))


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

ga_errs$err_type %>% unique
colnames(ga_errs)
understand_ga <- ga_errs %>%
  select(code.nodes, code.linker1, code.linker2, cat, err_type) %>% 
  #left_join(tob_L_categories, by="code.linker") %>% 
  mutate_all(funs(factor))  # set as factor so we can easily run summaries


# Honestly, after the analysis below, I think a bar graph might be the clearest way to make our point,
# followed by "see text/SI for discussion of classes of error". (specific examples of visualized CIFs)
# We just need to reorganize the bars, add percentages, and possibly color the bars by success vs.
# node incompatibility and other common sources of error, etc.
  
# Now let's dig into this data more closely
# For the figure in the paper, recall that my goal is to clearly show why these errors occur.
# We don't necessarily need all the detail, so an "Other" column is perfectly acceptable.
RUN_ERROR_FLOWS <- TRUE  # warning: disabling will also disable attribution notation
understand_tobacco <- understand_tobacco %>% mutate(err_cause = NA)
understand_ga <- understand_ga %>% mutate(err_cause = NA)
if (RUN_ERROR_FLOWS) {
  # There are two ways we can think about diagnosing the errors in the code
  # 1. Going from errors to their cause (this block)
  # 2. Going from patterns in the heatmap (e.g. certain linkers) back to types of error
  # Both approaches could be informative, but let's start with #1 here.
  
  # By category:
  # Trivial example summary for tpt vs. stp:
  understand_tobacco %>% filter(err_type=="Definition mismatch") %>% summary()
  understand_tobacco <- understand_tobacco %>% 
    mutate(err_cause = ifelse(err_type=="Definition mismatch" & code.topology=="tpt", "stp-tpt mixup", err_cause))
  # Most of the topological errors are sym_4_mc_1 (tetrahedral zinc).  Same with crashes
  understand_tobacco %>% filter(err_type=="topology") %>% summary()
  understand_tobacco %>% filter(err_type=="Crash") %>% summary()
  understand_tobacco <- understand_tobacco %>% 
    mutate(err_cause = ifelse(
      err_type %in% c("topology", "Crash") & mc1=="sym_4_mc_1",
      "sym_4_mc_1", err_cause
      ))
  # nodes errors are mostly sym_8_mc_7.  When including results from mc1/2, it's 273/282
  understand_tobacco %>% filter(err_type=="node_single_bonds") %>% summary()
  understand_tobacco <- understand_tobacco %>% 
    mutate(err_cause = ifelse(
      err_type=="node_single_bonds" & (mc1=="sym_8_mc_7" | mc2=="sym_8_mc_7"),
      "sym_8_mc_7", err_cause
      ))
  # "Many errors" is largely L_29 (strangely a benzene with double bonds coming off it)
  # Ah, and again sym_8_mc_7 problems propagating back to other sources of error 
  understand_tobacco %>% filter(err_type=="Many errors") %>% summary()
  tobacco_errs %>% filter(err_type=="Many errors" & code.linker=="L_29")
  understand_tobacco <- understand_tobacco %>% 
    mutate(err_cause = ifelse(
      err_type %in% c("Many errors", "Two errors") & code.linker=="L_29",
      "L_29", err_cause
    )) %>% 
    mutate(err_cause = ifelse(
      err_type %in% c("Many errors", "Two errors") & (mc1=="sym_8_mc_7" | mc2=="sym_8_mc_7"),
      "sym_8_mc_7", err_cause
    ))
  # TODO: Based on "two errors" and the common problems with sym_8_mc_7, let's make a new class of error
  # aggregating node and linker BO problems.
  # But again, half of the errors are sym_8_mc_7.
  understand_tobacco %>% filter(err_type=="Two errors") %>% summary
  tobacco_errs %>% filter(err_type=="Two errors")
  # Looking into the linkers, nothing particularly stands out from SB's
  understand_tobacco %>% filter(err_type=="linker_single_bonds") %>% summary
  # But 564/841 of the bond order problems are these 4-N rings
  understand_tobacco %>% filter(err_type=="linker_bond_orders") %>% summary
  understand_tobacco <- understand_tobacco %>% 
    mutate(err_cause = ifelse(
      err_type=="linker_bond_orders" & as.logical(L_4n_ring),
      "4-N ring", err_cause
    ))
  # Interestingly, all of the "formula" errors are L_8,9,10, which have big rings next to a phenyl.
  understand_tobacco %>% filter(err_type=="formula") %>% summary
  understand_tobacco <- understand_tobacco %>% 
    mutate(err_cause = ifelse(
      err_type=="formula" & code.linker %in% c("L_8", "L_9", "L_10"),
      "Overlapping rings", err_cause
    ))
  
  ### Now repeating for the GA MOFs ###
  
  # From the bar graph, crashes are the most important to address.
  # Interestingly, 101/129 are short L_0,1,2 secondary linkers
  understand_ga %>% filter(err_type=="Crash") %>% summary()
  short_ga <- c(0,1,2)
  understand_ga <- understand_ga %>% 
    mutate(err_cause = ifelse(
      err_type == "Crash" & (code.linker1 %in% short_ga | code.linker2 %in% short_ga),
      "Short linker (0-3)", err_cause
    ))
  understand_ga %>% filter(err_type=="Crash") %>% filter(!(code.linker1 %in% short_ga) & !(code.linker2 %in% short_ga)) %>% summary()
  # All of the definition mismatches are from Zr MOFs
  understand_ga %>% filter(err_type=="Definition mismatch") %>% summary()
  understand_ga <- understand_ga %>% 
    mutate(err_cause = ifelse(
      err_type == "Definition mismatch" & code.nodes == 4,
      "Zr coordination", err_cause
    ))
  # Many of the "two error" cases are the definition mismatch plus something else
  understand_ga %>% filter(err_type=="Two errors") %>% summary()
  ga_errs %>% filter(err_type=="Two errors")
  # Most of the single bonds errors are the gigantic L_4 linker
  understand_ga %>% filter(err_type=="linker_single_bonds") %>% summary()
  understand_ga <- understand_ga %>% 
    mutate(err_cause = ifelse(
      err_type == "linker_single_bonds" & (code.linker1==4 | code.linker2==4),
      "Bulky L_4", err_cause
    ))
  # 32/51 of the linker_bond_orders errors are the 4N linker L_27
  understand_ga %>% filter(err_type=="linker_bond_orders") %>% summary()
  understand_ga <- understand_ga %>% 
    mutate(err_cause = ifelse(
      err_type == "linker_bond_orders" & (code.linker1==27 | code.linker2==27),
      "4-N ring", err_cause
    ))
  # Formula errors?  Weirdly enough, the first linker is always either L_6 or L_27
  understand_ga %>% filter(err_type=="formula") %>% summary()
  #library(readr)
  #ga_errs %>% filter(err_type=="formula") %>% select(from_mofid) %>% write_delim("test.smi", col_names=FALSE, delim="|")
  # Even weirder, there is not a single correct 6/14 combo
  understand_ga <- understand_ga %>% 
    mutate(err_cause = ifelse(
      err_type == "formula" & (code.linker1==6 | code.linker2==14),
      "L_6/L_14 combo", err_cause
    )) %>% 
    mutate(err_cause = ifelse(
      err_type == "formula" & (code.linker1==27),
      "4-N ring", err_cause
    ))
  # Nonplanar carboxylates: all of the planar 4-c nodular "linkers" have this problem
  understand_ga %>% filter(err_type=="nonplanar_carboxylate") %>% summary()
  understand_ga <- understand_ga %>% 
    mutate(err_cause = ifelse(
      err_type == "nonplanar_carboxylate" & (code.linker1 %in% c(38,39)),
      "Planar 4-c 'linker'", err_cause
    ))
  # Overall, only 389 non-successes in this DB, including the definition mismatches
  understand_ga %>% filter(err_type!="Success") %>% nrow
  
  
  
  understand_tobacco <- understand_tobacco %>% mutate(err_cause = ifelse(is.na(err_cause), "Other", err_cause))
  understand_ga <- understand_ga %>% mutate(err_cause = ifelse(is.na(err_cause), "Other", err_cause))
  
}  # endif(RUN_ERROR_FLOWS)


### PLOTS FOR MAIN TEXT ###

error_code_translator <- tribble(
  ~err_type, ~err_plot,
  "Crash", "Program error",
  "linker_bond_orders", "Wrong bond orders",
  "linker_single_bonds", "Linker bonding",
  "node_single_bonds", "Node bonding",
  "Two errors", "Multiple errors",
  "Many errors", "Multiple errors",
  "topology", "Different topology",
  "formula", "Missing/extra atoms",
  "no_mof", "No MOF found",
  "nonplanar_carboxylate", "Nonplanar carboxylate",
  "Definition mismatch", "Misleading definition"
)

# Make error colors consistent between the two plots
all_err_causes <- 
  c(unique(understand_ga$err_cause), unique(understand_tobacco$err_cause)) %>% 
  unique() %>% 
  factor() %>% 
  levels()
# WARNING: the labels will break if more levels are added/removed, so let's add an assertion right here:
assert_that(length(all_err_causes) == 12)
err_colors_mapped <- brewer.pal(12, "Set3")
names(err_colors_mapped) <- all_err_causes
understand_ga <- understand_ga %>%
  mutate(err_cause = factor(err_cause, levels = all_err_causes)) %>% 
  mutate(err_cause = fct_drop(fct_infreq(err_cause)))
understand_tobacco <- understand_tobacco %>%
  mutate(err_cause = factor(err_cause, levels = all_err_causes)) %>% 
  mutate(err_cause = fct_drop(fct_infreq(err_cause)))
tobacco_err_names <- levels(understand_tobacco$err_cause)
tobacco_err_colors <- err_colors_mapped[tobacco_err_names]
ga_err_names <- levels(understand_ga$err_cause)
ga_err_colors <- err_colors_mapped[ga_err_names]


# Top bar for the overall success rate
plot_successes <- function(x, db_name) {
  num_mofs <- nrow(x)
  x %>% 
    mutate(success = ifelse(err_type=="Success", "Match", "Mismatch"), fake="MOFs") %>% 
    mutate(success = factor(success, levels=c("Mismatch", "Match"))) %>% 
    group_by(success, fake) %>% 
    summarize(success_rate = n()) %>% 
    mutate(match_label = paste0(success_rate, "\n", success)) %>% 
    ggplot(aes(x=fake, y=success_rate, fill=success, label=match_label)) +
    geom_col() +
    coord_flip() +
    theme_nothing() +
    geom_text(position=position_stack(vjust=0.5)) +  # thanks to https://stackoverflow.com/questions/6644997/showing-data-values-on-stacked-bar-chart-in-ggplot2
    ggtitle(paste(nrow(x), db_name, "in test set")) +
    theme(plot.title=element_text())
}

plot_curve_brace <- function(direction=1) {
  # Set direction of the sigmoid to +1 or -1 for upside down
  seq(0.0, 1.0, 0.001) %>% 
    data.frame(x=., y=sapply(., function(x) {direction*log((x)/(1-x))})) %>% 
    ggplot(aes(x, y)) +
    geom_point(color="#f8756d") +
    theme_nothing()
}

p_errors_tob <-
  understand_tobacco %>% 
  filter(err_type != "Success") %>% 
  left_join(error_code_translator, by="err_type") %>% 
  ggplot(aes(fct_rev(fct_infreq(err_plot)), fill=err_cause)) +
  geom_bar() +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values=tobacco_err_colors, breaks=tobacco_err_names) +
  coord_flip() +
  labs(x = NULL, y = NULL, fill = "Likely cause") +
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"))

p_errors_ga <-
  understand_ga %>% 
  filter(err_type != "Success") %>% 
  left_join(error_code_translator, by="err_type") %>% 
  ggplot(aes(fct_rev(fct_infreq(err_plot)), fill=err_cause)) +
  geom_bar() +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values=ga_err_colors, breaks=ga_err_names) +
  coord_flip() +
  labs(x = NULL, y = NULL, fill = "Likely cause") +
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"))

default_subfig_labels <- paste0("(", letters, ")")
plot_grid(
  # Adjust these brace locations as necessary for aesthetics:
  ggdraw() +
    draw_plot(plot_successes(understand_ga, "GA hMOFs"), y=0.8, height=0.2) +
    draw_plot(plot_curve_brace(1), y=0.74, height=0.07, x=0.0, width=0.90) +
    draw_plot(plot_curve_brace(-1), y=0.74, height=0.07, x=0.95, width=0.05) +
    draw_plot(p_errors_ga, y=0.0, height=0.75),
  ggdraw() +
    draw_plot(plot_successes(understand_tobacco, "ToBaCCo MOFs"), y=0.8, height=0.2) +
    draw_plot(plot_curve_brace(1), y=0.74, height=0.07, x=0.0, width=0.80) +
    draw_plot(plot_curve_brace(-1), y=0.74, height=0.07, x=0.90, width=0.10) +
    draw_plot(p_errors_tob, y=0.0, height=0.75),
  nrow = 2,
  labels = default_subfig_labels, label_size=14, hjust=0, vjust=1
  ) %>% 
  cowplot::save_plot(
    "Analysis/Figures/roundtrip.png", .,
    base_aspect_ratio=1.5, nrow = 2
  )

