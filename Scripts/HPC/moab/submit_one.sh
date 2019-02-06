#!/bin/bash
# Call as Scripts/HPC/submit_one.job path_to_cifs short_friendly_job_name

dir=$1
echo "Submitting $dir"
msub Scripts/HPC/child_mofid_moab.job -N id-$2 -v CIFDIR=$dir,CIFOUTLABEL=$2,OUTPUT_DIR="Output_$2",MOFSTDOUT="out_$2.smi" -o pbs_out_$2.txt -e err_$2.txt

