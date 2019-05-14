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

COPY_MOFID="${OUTPUT_DIR}/folder_mofid.smi"
COPY_MOFKEY="${OUTPUT_DIR}/folder_mofkey.tsv"

SRC_LINKER_STATS="${OUTPUT_DIR}/MetalOxo/linker_stats.txt"
COPY_LINKER_STATS="${OUTPUT_DIR}/folder_linker_stats.tsv"

rm -f "${COPY_MOFID}"
echo -e "filename\tmofkey" > "${COPY_MOFKEY}"
echo -e "filename\tinchikey\tconnections_metaloxo_net\tuc_count\tinchi\ttruncated_inchikey\tsmiles\tskeleton" > "${COPY_LINKER_STATS}"

echo "Analyzing ${CIF_DIR} with MOFid commit:" 1>&2
git rev-parse --verify HEAD 1>&2
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
	
	# Also parse the linker stats (deleting blank lines, first)
	sed -e '/^$/d' "${SRC_LINKER_STATS}" | sed -e 's/^/'"$(basename "$i")"'\t/' >> "${COPY_LINKER_STATS}"
	rm -f "${SRC_LINKER_STATS}"
done

echo "----------------------------" 1>&2
echo "Time:" 1>&2
date 1>&2
echo "----------------------------" 1>&2
