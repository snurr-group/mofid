#!/bin/bash
# Gathers job outputs on NUIT quest from run_submit_all.sh

SUMMARY_DIR="Summary/"
ORIG_OUTPUTS="$SUMMARY_DIR/Data/"
JOB_SCHEDULER="moab"


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


mv Output_*/ "$ORIG_OUTPUTS"
mv err_*.txt "$ORIG_OUTPUTS"
mv pbs_out_*.txt "$ORIG_OUTPUTS"
mv out_*.smi "$ORIG_OUTPUTS"
mv jobid_* "$ORIG_OUTPUTS"


echo "Submitting jobs to validate MOFid's ToBaCCo MOFs and GA hMOFs against their recipes..."
if [ "$JOB_SCHEDULER" == "moab" ]
then
	msub Scripts/HPC/moab/child_validation.job -N validate-ga -v BASE_SMILES=ga,SUMMARY_DIR=${SUMMARY_DIR} -o ${ORIG_OUTPUTS}/pbs_out_ga_validation.txt -e ${ORIG_OUTPUTS}/err_ga_validation.txt
	sleep 5
	msub Scripts/HPC/moab/child_validation.job -N validate-tob -v BASE_SMILES=tob,SUMMARY_DIR=${SUMMARY_DIR} -o ${ORIG_OUTPUTS}/pbs_out_tob_validation.txt -e ${ORIG_OUTPUTS}/err_tob_validation.txt
	sleep 5
elif [ "$JOB_SCHEDULER" == "slurm" ]
then
	echo "slurm not yet supported.  See other scripts for an example for how to adapt moab" 1>&2
	exit 2
else
	echo "Unknown JOB_SCHEDULER.  Cannot run validation on the output.smi files" 1>&2
	exit 2
fi

