#!/bin/bash
# Call as Scripts/HPC/moab/submit_parallel.job path_to_cifs_parent_dir_of_cifx short_friendly_job_name

for dir in cif{1..5}
do
	echo "Submitting $dir"
	msub Scripts/HPC/moab/child_mofid_moab.job -N id-$2-$dir -v CIFDIR=$1/$dir,CIFOUTLABEL=$2,OUTPUT_DIR="Output_$2_$dir",MOFSTDOUT="out_$2_$dir.smi" -o pbs_out_$2_$dir.txt -e err_$2_$dir.txt
	sleep 10
done
