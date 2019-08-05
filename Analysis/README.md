# Analysis for MOFid paper

Analysis scripts to generate figures and stats for the manuscript


## Steps to reproduce figures

1. Copy relevant MOF databases (GA hMOFs, ToBaCCo, and CoRE MOF 2019) to the Data/ directory.
2. Generate MOFid's and validation data on an HPC cluster, e.g. [Northwestern quest](https://www.it.northwestern.edu/research/user-services/quest/acknowledgment.html).
	1. `Scripts/HPC/run_submit_all.sh` and wait for jobs to finish.
	2. `Scripts/HPC/run_post_submit_all.sh` to collect results and run the validator.
3. Extract the carbon-containing CoRE MOFs by running `python Analysis/extract_carbon_mofs.py`
4. Convert `Summary` results to the corresponding tsv files: `python Python/convert_smi_to_tables.py Summary/core.smi && Scripts/combine_smi_tables.sh`
	1. The bash scripts require [sqlite3](https://www.sqlite.org/index.html) to be installed and in your path.
5. Open the `mofid.Rproj` file in RStudio, then run analysis scripts to generate validation figures.
	1. `Analysis/plot_round_trip_errors.R` parses the validator .json output to find a percent success/match rate for the MOFid code vs. structure recipe, as well as figures showing the composition.
6. Copy files as necessary to the Figures subdirectories from Summary/
7. `cd Figures/Polymorphs && ../../find_polymorphs.sh combined_table.tsv out_with_names.tsv out_summary.tsv`
8. `cd Figures/Duplicates; for db in core tob ga; do ../../find_duplicates.sh duplicates ${db}_mofkey.tsv duplicates_${db}_names.tsv duplicates_${db}_summary.tsv duplicates_${db}_all_families.tsv; done;`
9. `cd Figures/Duplicates && ../../find_duplicates.sh overlap core_mofkey.tsv ga_mofkey.tsv overlap_with_ga_names.tsv overlap_with_ga_summary.tsv && ../../find_duplicates.sh overlap core_mofkey.tsv tob_mofkey.tsv overlap_with_tob_names.tsv overlap_with_tob_summary.tsv && ../../find_duplicates.sh overlap tob_mofkey.tsv ga_mofkey.tsv overlap_tob_and_ga_names.tsv overlap_tob_and_ga_summary.tsv && ../../find_duplicates.sh overlap common_mofkeys.tsv core_mofkey.tsv overlap_common_names.tsv overlap_common_summary.tsv`
10. Run `plot_venn_diagram.R`, which uses the results saved in `Figures/Duplicates`
11. Copy Data/Output_core_cif*/folder_smiles_parts*.tsv to FolderData (according to the template), then run the common node/linker analysis via: `cd Figures/Database-Stats/ && cp ../../../Summary/OrigCore/* CopiedOrigCore && cp ../../../Summary/core* . && cd FolderData && ./combine_folders.sh && cd .. && mkdir -p Output && ./run_stats.sh unfiltered_folder_smiles_parts.tsv CopiedOrigCore/core_linker_stats.tsv core_extracted_files_with_carbon.txt Output/parts_summary.tsv Output/parts_linker25.smi Output/parts_nodes.tex Output/stats_summary.tsv Output/stats_linker20.smi`
12. See also the notebooks in `Figures/` subdirectories


## Other comments

Analysis scripts will live here.  More general-purpose utilities are in the /bin, /Python, and /Scripts directories.

