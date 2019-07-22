#!/bin/bash
# Runs the Python MOFid analysis on a directory of CIFs
#
# Outputs the MOFids in the .smi SMILES format to stdout and any logging
# info, such as warnings from sbu.exe or the current progress, to stderr.
#
# The intermediate OUTPUT_DIR is temporary and overwritten by each MOF,
# but it is included as a parameter to allow multiple copies of run_folder.sh
# to run concurrently, such as parallel jobs on an HPC cluster.  Otherwise,
# multiple processes would be writing to the same intermediate topology folder
# and possibly cause inconsistencies.

# Parse input arguments
CIF_DIR="$1"
OUTPUT_DIR="$2"
if [ $# -ne 2 ]
then
	echo "ERROR: Incorrect number of arguments" 1>&2 && exit
fi

OUTPUT_COPY="${OUTPUT_DIR}/results_part.txt"
mkdir -p "${OUTPUT_DIR}"
PYTHON_MOFID="${OUTPUT_DIR}/python_mofid.txt"
PYTHON_MOFKEY="${OUTPUT_DIR}/python_mofkey.txt"
PYTHON_SMILES_PARTS="${OUTPUT_DIR}/python_smiles_parts.txt"
PYTHON_MOLEC_FORMULA="${OUTPUT_DIR}/python_molec_formula.txt"

COPY_MOFID="${OUTPUT_DIR}/folder_mofid.smi"
COPY_MOFKEY="${OUTPUT_DIR}/folder_mofkey.tsv"
COPY_SMILES_PARTS="${OUTPUT_DIR}/folder_smiles_parts.tsv"
COPY_MOLEC_FORMULA="${OUTPUT_DIR}/folder_molec_formula.tsv"

SRC_LINKER_STATS="${OUTPUT_DIR}/MetalOxo/linker_stats.txt"
COPY_LINKER_STATS="${OUTPUT_DIR}/folder_linker_stats.tsv"

rm -f "${COPY_MOFID}"
echo -e "filename\tmofkey" > "${COPY_MOFKEY}"
echo -e "filename\tinchikey\tconnections_metaloxo_net\tuc_count\tinchi\ttruncated_inchikey\tsmiles\tskeleton" > "${COPY_LINKER_STATS}"
echo -e "filename\tpart\tsmiles" > "${COPY_SMILES_PARTS}"
echo -e "filename\tdelim_formula" > "${COPY_MOLEC_FORMULA}"

# Save current git commit to improve traceability
echo "Analyzing ${CIF_DIR} with MOFid commit:" 1>&2
# git rev-parse --verify HEAD 1>&2
# Replace explicit dependence on git with parsing the relevant .git files
# Thanks to Fordi: https://stackoverflow.com/questions/949314/how-to-retrieve-the-hash-for-the-current-commit-in-git/33133769#33133769 
HASH="ref: HEAD"; while [[ $HASH == ref\:* ]]; do HASH="$(cat ".git/$(echo $HASH | cut -d \  -f 2)")"; done;
echo $HASH 1>&2

echo "Writing results to intermediate directory: ${OUTPUT_DIR}" 1>&2
echo "----------------------------" 1>&2
echo "Time:" 1>&2
date 1>&2
echo "----------------------------" 1>&2

for i in "${CIF_DIR}/"*.cif
do
	echo "Beginning $i" 1>&2
	echo "----------------------------" 1>&2
	python Python/run_mofid.py "$i" "${OUTPUT_DIR}" | tee -a "${OUTPUT_COPY}"
	
	# Make an extra copy of the MOFid and MOFkey as written to files
	cat "${PYTHON_MOFID}" >> "${COPY_MOFID}"
	echo -e "$(basename "$i")\t$(tr -d '\n' < "${PYTHON_MOFKEY}")" >> "${COPY_MOFKEY}"
	rm -f "${PYTHON_MOFID}" "${PYTHON_MOFKEY}"
	
	# Also parse the linker stats (deleting blank lines, first) and SMILES parts
	sed -e '/^$/d' "${SRC_LINKER_STATS}" | sed -e 's/^/'"$(basename "$i")"'\t/' >> "${COPY_LINKER_STATS}"
	sed -e '/^$/d' "${PYTHON_SMILES_PARTS}" | sed -e 's/^/'"$(basename "$i")"'\t/' >> "${COPY_SMILES_PARTS}"
	sed -e '/^$/d' "${PYTHON_MOLEC_FORMULA}" | sed -e 's/^/'"$(basename "$i")"'\t/' >> "${COPY_MOLEC_FORMULA}"
	rm -f "${SRC_LINKER_STATS}" "${PYTHON_SMILES_PARTS}" "${PYTHON_MOLEC_FORMULA}"
done

echo "----------------------------" 1>&2
echo "Time:" 1>&2
date 1>&2
echo "----------------------------" 1>&2
