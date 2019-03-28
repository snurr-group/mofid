# Analysis for MOFid paper

Analysis scripts to generate figures and stats for the manuscript


## Steps to reproduce figures

1. Copy relevant MOF databases (GA hMOFs, ToBaCCo, and CoRE MOF 2019) to the Data/ directory.
2. Generate MOFid's and validation data on an HPC cluster, e.g. [Northwestern quest](https://www.it.northwestern.edu/research/user-services/quest/acknowledgment.html).
	1. `Scripts/HPC/run_submit_all.sh` and wait for jobs to finish.
	2. `Scripts/HPC/run_post_submit_all.sh` to collect results and run the validator.
3. TODO: probably something with `Python/convert_smi_to_tables.py` and `Scripts/combine_smi_tables.sh`
	1. The bash scripts require [sqlite3](https://www.sqlite.org/index.html) to be installed and in your path.
	2. TODO: Still need to implement find_duplicate_mofids.sh and similar for database_overlaps, but consider whether we do this analysis via sqlite3, how to encode group numbers/parent database (like Barthel et al.?), etc.
4. Open the `mofid.Rproj` file in RStudio, then run analysis scripts to generate validation figures.
	1. `Analysis/plot_round_trip_errors.R` parses the validator .json output to find a percent success/match rate for the MOFid code vs. structure recipe, as well as figures showing the composition.
5. TODO: other scripts in this directory for analysis


## Other comments

Analysis scripts will live here.  More general-purpose utilities are in the /bin, /Python, and /Scripts directories.

