# Analysis for MOFid paper

Analysis scripts to generate figures and stats for the manuscript


## Steps to reproduce figures

1. Copy relevant MOF databases (GA hMOFs, ToBaCCo, and CoRE MOF 2019) to the Data/ directory.
2. Generate MOFid's and validation data on an HPC cluster, e.g. [Northwestern quest](https://www.it.northwestern.edu/research/user-services/quest/acknowledgment.html).
	1. `Scripts/HPC/run_submit_all.sh` and wait for jobs to finish.
	2. `Scripts/HPC/run_post_submit_all.sh` to collect results and run the validator.
3. TODO: probably something with `Python/convert_smi_to_tables.py`
4. TODO: scripts in this directory for analysis


## Other comments

Analysis scripts will live here.  More general-purpose utilities are in the /bin, /Python, and /Scripts directories.

