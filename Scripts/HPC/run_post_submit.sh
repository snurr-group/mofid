#!/bin/bash
# Gathers job outputs on NUIT quest from run_submit_all.sh

SUMMARY_DIR="Summary/"
ORIG_OUTPUTS="$SUMMARY_DIR/Data/"
JOB_SCHEDULER="slurm"


if [ -d "$SUMMARY_DIR" ]
then
	echo "Output directory $SUMMARY_DIR already exists!  Exiting without making changes" 1>&2
	exit 1
fi

mkdir -p "$SUMMARY_DIR"
mkdir -p "$ORIG_OUTPUTS"

cat err_tob_cif{1..5}.txt > "$SUMMARY_DIR/err_tob.txt"
cat err_core_cif{1..5}.txt > "$SUMMARY_DIR/err_core.txt"
cp err_ga.txt "$SUMMARY_DIR/err_ga.txt"

cat out_tob_cif{1..5}.smi > "$SUMMARY_DIR/tob.smi"
cat out_core_cif{1..5}.smi > "$SUMMARY_DIR/core.smi"
cp out_ga.smi "$SUMMARY_DIR/ga.smi"


# Also copy MOFid, MOFkey, and linker stats
cat Output_core_cif{1..5}/folder_mofid.smi > "$SUMMARY_DIR/core_mofid.smi"
cat Output_tob_cif{1..5}/folder_mofid.smi > "$SUMMARY_DIR/tob_mofid.smi"
cp Output_ga/folder_mofid.smi "$SUMMARY_DIR/ga_mofid.smi"
copy_folder_tsv() {
	# Copies files named folder_SUFFIX.tsv for core/tob/ga
	# Use the `tail -q` flag to avoid extra filename headers in the output
	suffix_name=$1
	head -n 1 Output_core_cif1/folder_${suffix_name}.tsv > "$SUMMARY_DIR/core_${suffix_name}.tsv"
	tail -q -n +2 Output_core_cif{1..5}/folder_${suffix_name}.tsv >> "$SUMMARY_DIR/core_${suffix_name}.tsv"
	head -n 1 Output_tob_cif1/folder_${suffix_name}.tsv > "$SUMMARY_DIR/tob_${suffix_name}.tsv"
	tail -q -n +2 Output_tob_cif{1..5}/folder_${suffix_name}.tsv >> "$SUMMARY_DIR/tob_${suffix_name}.tsv"
	cp Output_ga/folder_${suffix_name}.tsv "$SUMMARY_DIR/ga_${suffix_name}.tsv"
}
copy_folder_tsv "mofkey"
copy_folder_tsv "linker_stats"
copy_folder_tsv "smiles_parts"
copy_folder_tsv "molec_formula"

mv Output_*/ "$ORIG_OUTPUTS"
mv err_*.txt "$ORIG_OUTPUTS"
mv pbs_out_*.txt "$ORIG_OUTPUTS"
mv out_*.smi "$ORIG_OUTPUTS"
mv jobid* "$ORIG_OUTPUTS"
mv slurm-*.out "$ORIG_OUTPUTS"
mkdir -p "$ORIG_OUTPUTS/Cores"
mv core.* "$ORIG_OUTPUTS/Cores/"


echo "Submitting jobs to validate MOFid's ToBaCCo MOFs and GA hMOFs against their recipes..."
if [ "$JOB_SCHEDULER" == "moab" ]
then
	msub Scripts/HPC/moab/child_validation.job -N validate-ga -v BASE_SMILES=ga,SUMMARY_DIR=${SUMMARY_DIR} -o ${ORIG_OUTPUTS}/pbs_out_ga_validation.txt -e ${ORIG_OUTPUTS}/err_ga_validation.txt
	sleep 5
	msub Scripts/HPC/moab/child_validation.job -N validate-tob -v BASE_SMILES=tob,SUMMARY_DIR=${SUMMARY_DIR} -o ${ORIG_OUTPUTS}/pbs_out_tob_validation.txt -e ${ORIG_OUTPUTS}/err_tob_validation.txt
	sleep 5
elif [ "$JOB_SCHEDULER" == "slurm" ]
then
	export BASE_SMILES=ga
	export SUMMARY_DIR=${SUMMARY_DIR}
	sbatch -J validate-ga -o ${ORIG_OUTPUTS}/pbs_out_ga_validation.txt -e ${ORIG_OUTPUTS}/err_ga_validation.txt Scripts/HPC/slurm/child_validation.job
	sleep 5
	export BASE_SMILES=tob
	sbatch -J validate-tob -o ${ORIG_OUTPUTS}/pbs_out_tob_validation.txt -e ${ORIG_OUTPUTS}/err_tob_validation.txt Scripts/HPC/slurm/child_validation.job
	sleep 5
else
	echo "Unknown JOB_SCHEDULER.  Cannot run validation on the output.smi files" 1>&2
	exit 2
fi

