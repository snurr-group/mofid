#!/bin/bash
# Submits all relevant jobs for NUIT quest HPC cluster.
# Generates output files to post-process with run_post_submit.sh
# to get MOFid's for CoRE MOF 2019-ASR, ToBaCCo, and a GA hMOF test subset

JOB_SCHEDULER="moab"

echo "Submitting ToBACCo MOFs"
Scripts/HPC/${JOB_SCHEDULER}/submit_parallel.sh Data/opt_tobacco/split-all/ tob
echo "Submitting CoRE MOF 2019"
Scripts/HPC/${JOB_SCHEDULER}/submit_parallel.sh Data/core2019 core
echo "Submitting GA hMOFs"
Scripts/HPC/${JOB_SCHEDULER}/submit_one.sh Data/CatParentGA/ ga
echo "Done submitting HPC jobs."

