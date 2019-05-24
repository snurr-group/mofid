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


# What are the identities of these MOFs?
# e.g. o_core_ga %>% anti_join(o_combined, by="identifier") %>% View
o_identities <- tribble(
  ~identifier, ~common_name, ~ref,
  # Overlap between all databases
  "Cu.QMKYBPDZANOJGF.MOFkey-v1.tbo", "Cu-BTC", "TODO",
  "Zn.KKEYFWRCBNTPAC.MOFkey-v1.pcu", "MOF-5", "TODO",
  "Zn.RXOHFPCZGPKIRD.MOFkey-v1.pcu", "IRMOF-8", "TODO",
  "Zr.KKEYFWRCBNTPAC.MOFkey-v1.fcu", "UiO-66", "TODO",
  "Zn.NEQFBGHQPUXOFH.MOFkey-v1.pcu", "IRMOF-10", "TODO",
  "Zn.ARBAMECCCQZCSR.MOFkey-v1.pcu", "IRMOF-61 (FAHPOV)", "10.1021/jp302356q (orig? 10.1039/C1CE06044A)",
  "Cu.NIJMZBWVSRWZFZ.MOFkey-v1.tbo", "TCM-8", "10.1039/C3CC49829H",
  "Cu.SATWKVZGMWCXOJ.MOFkey-v1.tbo", "DUT-34", "10.1002/chem.201101383",
  # CoRE - GA
  "V.KKEYFWRCBNTPAC.MOFkey-v1.rna", "MIL-47", "TODO",
  "Zn.MWVTWFVJZLCBMC.NEQFBGHQPUXOFH.MOFkey-v1.pcu", "BMOF-1-bpdc", "10.1021/ic202683s",
  "Zn.QMKYBPDZANOJGF.MOFkey-v1.tbo", "Zn-HKUST-1", "10.1039/C2CE26115D",
  "V.NEQFBGHQPUXOFH.MOFkey-v1.rna", "VO(BPDC)", "10.1021/ic301338a",
  "Zn.KKEYFWRCBNTPAC.NEQFBGHQPUXOFH.MOFkey-v1.pcu", "SUMOF-4", "10.1039/C2JM15933C",
  "V.RXOHFPCZGPKIRD.MOFkey-v1.rna", "COMOC‐3", "10.1002/ejic.201101099",
  # ToBaCCo - GA: structures without a corresponding CoRE MOF
  # CoRE - ToBaCCo
  "Cu.SATWKVZGMWCXOJ.MOFkey-v1.pto", "MOF-143", "10.1021/ic201376t",
  "Zr.NEQFBGHQPUXOFH.MOFkey-v1.fcu", "UiO-67", "orig, but not in the WIZMAV? 10.1021/ja8057953",
  "Zn.NWYGETXZXGDGKD.MOFkey-v1.pyr", "SNU-77", "10.1002/chem.201003376",
  "Zr.VZCYOOQTPOCHFL.MOFkey-v1.fcu", "MOF-801", "10.1021/ja500330a and 10.1016/j.micromeso.2011.12.010",
  "Co.SATWKVZGMWCXOJ.MOFkey-v1.the", NA, "Maybe? 10.1039/C3SC51379C",
  "Cu.NIJMZBWVSRWZFZ.MOFkey-v1.pto", "TCM-4", "10.1002/chem.201304856",
  "Zr.HVCDAMXLLUJLQZ.MOFkey-v1.csq", "NU-1000", "TODO",
  "Zr.MSFXUHUYNSYIDR.MOFkey-v1.spn", "PCN-777, 493-MOF-BA", "10.1002/anie.201409334 and 10.1039/C7CC00029D",
  "Co.KKEYFWRCBNTPAC.MOFkey-v1.bcu", NA, "Maybe a MOF? 10.1021/ic900816h",
  "Cu.MSFXUHUYNSYIDR.MOFkey-v1.tbo", "PCN-6'", "10.1021/ja067435s",
  "Cu.OMMYHUPHWRFAPM.MOFkey-v1.lvt", "ZJU-30", "10.1021/cg400164m",
  "Cu.QURGMSIQFRADOZ.MOFkey-v1.nbo", "MOF-505", "10.1002/anie.200462787",
  "Cu.STBMCGYZKQKNCV.MOFkey-v1.lvt", "UTSA-60", "10.1039/C4CC09999K",
  "Zn.YROTZTMCXKTYMW.MOFkey-v1.pcu", "IRMOF-62", "10.1021/jp302356q and 10.1039/C1CE05354J",
  "Zr.ARBAMECCCQZCSR.MOFkey-v1.fcu", "BUT-30", "10.1016/j.jssc.2014.07.001",
  "Zr.FZTIWOBQQYPTCJ.MOFkey-v1.fcu", "UiO-68", "10.1038/ncomms12610",
  "Zr.IIIWRSPHUBZZOB.MOFkey-v1.csq", "PCN-128", "10.1021/jacs.5b04695",
  "Zr.KKCMOJXVSJSNKW.MOFkey-v1.ftw", "NU-1100", "10.1002/chem.201402895",
  "Co.ABMFBCRYHDZLRD.MOFkey-v1.bcu", NA, NA,
  "Co.MSFXUHUYNSYIDR.MOFkey-v1.the", "PCN-9 (Co)", "10.1021/ja063538z and 10.1021/ic900475q",
  "Co.VSFXBCHNPQPWBX.MOFkey-v1.flu", "SNU-15", "10.1039/B900085B",
  "Cu.ABMFBCRYHDZLRD.MOFkey-v1.lvt", "(SUJNUH)", "10.1039/B917029D",
  "Cu.FUVQBDCRAGNDAO.MOFkey-v1.pts", "(OWIZAW)", "10.1021/ic1009169",
  "Cu.NWYGETXZXGDGKD.MOFkey-v1.pto", "DUT-64", "10.1039/C6CE01513A",
  "Cu.NWYGETXZXGDGKD.MOFkey-v1.tbo", "DUT-63", "10.1039/C6CE01513A",
  "Cu.OMMYHUPHWRFAPM.MOFkey-v1.pts", "UTSA-68", "10.1039/C5CC10598F",
  "Cu.ONMZLXHMUBUCKZ.MOFkey-v1.pts", "(SUKYON)", "10.1002/anie.200904983",
  "Cu.OYJYITDWCXDYQZ.MOFkey-v1.tbo", "ZJU-199", "10.1021/acs.cgd.6b01382",
  "Cu.PEQRGMPXYDIZSX.MOFkey-v1.tbo", "MOF-399", "10.1021/ic201376t",
  "Cu.SVQPJYAUEQAOEU.MOFkey-v1.pts", "IMP-9", "10.1021/cg1008768",
  "Cu.VEBUOOBGPZWCFE.MOFkey-v1.pto", "Cu−TCA", "10.1002/adfm.201102157",
  "Cu.VSFXBCHNPQPWBX.MOFkey-v1.pts", "(SUKYIH)", "10.1002/anie.200904983",
  "Cu.VTROUXWDAQESAI.MOFkey-v1.pto", "UTSA-28-Cu", "10.1021/ic401870e",
  "Cu.WCKARIXPDNMKGB.MOFkey-v1.pto", "MOF-388 (uncatenated)", "10.1021/ic201376t",
  "Mn.BGMFKFGEYWEPCF.MOFkey-v1.the", "(BIWSIK)", "10.1021/ic701917w",
  "Mn.GRYXMEJPMFYXQG.MOFkey-v1.the", "(JEWYAM)", "10.1021/ja0656853",
  "Mn.UBTDCKNMTXZSND.MOFkey-v1.flu", "IMP-16Mn", "10.1039/C4CE00486H",
  "Zn.KMOBUNLHXZTQGA.MOFkey-v1.pyr", "SNU-150", "10.1002/chem.201303086",
  "Zn.OTAJGWQCQIEFEV.MOFkey-v1.pcu", "IRMOF-14", "10.1126/science.1067208",  # both IRMOF-14 and LIHFAK were in the original hMOFs!!!  Why not in this run?
  "Zn.SBBQDUFLZGOASY.MOFkey-v1.pcu", "(LIHFAK)", "10.1021/ja0700395",
  "Zn.VEBUOOBGPZWCFE.MOFkey-v1.pyr", "(NAPFUG)", "10.1021/ja043756x",
  "Zr.HVWAJJAMCRZUIQ.MOFkey-v1.ftw", "(MUBZOA)", "10.1002/anie.201406501",
  "Zr.KVQMUHHSWICEIH.MOFkey-v1.fcu", "(XIVTED)", "10.1039/C3CC48275H",
  "Zr.OMMYHUPHWRFAPM.MOFkey-v1.flu", "UMCM-312", "10.1021/acs.cgd.6b00698",
  "Zr.PLGNSZYEWBJPKZ.MOFkey-v1.ftw", "NU-1102", "10.1021/ja512973b",
  "Zr.QMKYBPDZANOJGF.MOFkey-v1.spn", "MOF-808", "10.1021/ja500330a",
  "Zr.SVQPJYAUEQAOEU.MOFkey-v1.flu", "(XUBGUY)", "10.1039/C5DT00421G",
  "Zr.UEGVXEAWMLARMQ.MOFkey-v1.ftw", "NU-1103", "10.1021/ja512973b",
  "Zr.VSFXBCHNPQPWBX.MOFkey-v1.flu", "MOF-841", "10.1021/ja500330a",
  "Zr.WVBWNURTSYQJFE.MOFkey-v1.ftw", "NU-1004", "10.1021/ja512973b",
  "Zr.YROTZTMCXKTYMW.MOFkey-v1.fcu", "(UKIBIB)", "10.1002/chem.201505185"
)

# TODO: considering a table/figure drawing the chemical structures?
# Could have boxes for the different nodes and categorize them that way


