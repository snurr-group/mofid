#!/bin/bash
# Call as Scripts/HPC/slurm/submit_parallel.sh path_to_cifs_parent_dir_of_cifx short_friendly_job_name

for dir in cif{1..5}
do
	echo "Submitting $dir"
	export CIFDIR=$1/$dir
	export CIFOUTLABEL=$2
	export OUTPUT_DIR="Output_$2_$dir"
	export MOFSTDOUT="out_$2_$dir.smi" 
	sbatch -J id-$2-$dir -o pbs_out_$2_$dir.txt -e err_$2_$dir.txt Scripts/HPC/slurm/child_mofid_slurm.job
	sleep 10
done
